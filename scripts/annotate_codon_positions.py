#!/usr/bin/env python3
"""
annotate_codon_positions.py

For each gene in a list:
  1. Loads the corresponding FASTA alignment
  2. Extracts the gene symbol from the Gene_Name: field in the alignment headers
  3. Fetches the canonical CDS from Ensembl Plants by gene symbol
  4. Finds the A. thaliana sequence (Species:Arabidopsis_thaliana) in the alignment
  5. Locates the CDS within the ungapped A. thaliana sequence
  6. Threads codon position annotations (1/2/3) back through the gapped alignment
  7. Writes a TSV annotation file per gene

Usage:
    python annotate_codon_positions.py \
        --genes gene_list.txt \
        --aln_dir /path/to/alignments \
        --out_dir /path/to/output \
        [--aln_suffix .dna.aln.fasta] \
        [--ath_species Arabidopsis_thaliana] \
        [--delay 0.3]

gene_list.txt: one numeric gene/family ID per line matching the alignment filename prefix
Alignment files must be named {id}{aln_suffix}, e.g. 4893.dna.aln.fasta
Headers must contain Gene_Name:<symbol> and Species:<name> fields.
"""

import argparse
from datetime import UTC, datetime
import logging
import re
import time
from pathlib import Path

import requests
from Bio import AlignIO

# ── Logging ──────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Ensembl Plants API ────────────────────────────────────────────────────────

ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}


def fetch_cds_ensembl(gene_symbol: str, session: requests.Session) -> tuple[str, str] | tuple[None, None]:
    """
    Fetch the canonical CDS for an A. thaliana gene.
    Uses /lookup/id/ for TAIR IDs, /lookup/symbol/ otherwise.
    Returns (CDS, transcript_id), or (None, None) on failure.
    """
    # TAIR IDs (e.g. AT1G49980) use the ID endpoint; gene symbols use the symbol endpoint
    if re.match(r'^AT[1-5MC]G\d+$', gene_symbol, re.IGNORECASE):
        url = f"{ENSEMBL_REST}/lookup/id/{gene_symbol}"
        params = {"expand": 1}
    else:
        url = f"{ENSEMBL_REST}/lookup/symbol/arabidopsis_thaliana/{gene_symbol}"
        params = {"expand": 1}

    try:
        r = session.get(url, headers=HEADERS, params=params, timeout=30)
        r.raise_for_status()
        gene_data = r.json()
    except Exception as e:
        log.warning(f"{gene_symbol}: gene lookup failed — {e}")
        return None, None

    # Find canonical transcript
    transcripts = gene_data.get("Transcript", [])
    if not transcripts:
        log.warning(f"{gene_symbol}: no transcripts found in Ensembl")
        return None, None

    # Prefer the transcript flagged is_canonical, else take the longest CDS
    canonical = next((t for t in transcripts if t.get("is_canonical") == 1), None)
    if canonical is None:
        canonical = max(
            transcripts,
            key=lambda t: t.get("Translation", {}).get("length", 0) if t.get("Translation") else 0,
        )

    transcript_id = canonical["id"]

    # Step 2: fetch CDS sequence for that transcript
    url = f"{ENSEMBL_REST}/sequence/id/{transcript_id}"
    try:
        r = session.get(
            url,
            headers=HEADERS,
            params={"type": "cds", "multiple_sequences": 0},
            timeout=30,
        )
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        log.warning(f"{gene_symbol} (transcript {transcript_id}): CDS fetch failed — {e}")
        return None, None

    if isinstance(data, list):
        data = data[0]
    cds = data.get("seq", "").upper()
    if not cds:
        log.warning(f"{gene_symbol}: empty CDS returned")
        return None, None

    log.info(f"{gene_symbol}: CDS length {len(cds)} bp (transcript {transcript_id})")
    return cds, transcript_id


# ── Alignment helpers ─────────────────────────────────────────────────────────


def extract_gene_name(alignment) -> str:
    """Parse Gene_Name: value from the first alignment record's description."""
    for rec in alignment:
        m = re.search(r'Gene_Name:(\S+)', rec.description)
        if m:
            return m.group(1)
    raise ValueError("Gene_Name field not found in alignment headers")


def find_ath_record(alignment, ath_species: str):
    """
    Return the A. thaliana record by matching ath_species in the full description.
    Raises ValueError if not found.
    """
    for rec in alignment:
        if ath_species in rec.description:
            return rec
    sample = [r.description[:80] for r in list(alignment)[:5]]
    raise ValueError(
        f"Cannot find sequence with '{ath_species}' in alignment. "
        f"Sample headers: {sample}"
    )


def ungap(seq_str: str) -> str:
    return seq_str.replace("-", "").replace(".", "").upper()


def find_cds_in_seq(seq: str, cds: str) -> tuple[str, int, int] | tuple[None, None, None]:
    """
    Locate the relationship between the alignment sequence and the CDS.

    Returns (mode, start, end) where:
      mode='flanked' — CDS is a substring of seq (alignment has UTR/flanking)
        start/end = CDS boundaries in ungapped seq coordinates
      mode='trimmed' — seq is a substring of CDS (alignment is a CDS subset)
        start/end = where seq begins in the CDS (used to compute codon offset)
      mode='full'    — lengths match within 3 bp; treat whole alignment as CDS
        start=0, end=len(seq)

    Returns (None, None, None) if no match found.
    """
    # Case 1: alignment contains the full CDS (has flanking UTR/sequence)
    for query in (cds, cds[:-3] if len(cds) > 3 else None):
        if query is None:
            continue
        idx = seq.find(query)
        if idx != -1:
            if query != cds:
                log.debug("CDS matched after trimming stop codon")
            return "flanked", idx, idx + len(query)

    # Case 2: alignment is a trimmed subset of the CDS
    idx = cds.find(seq)
    if idx != -1:
        return "trimmed", idx, idx + len(seq)

    # Case 3: lengths within 3 codons — likely same gene, different version
    if abs(len(seq) - len(cds)) <= 9:
        return "near_full", 0, len(seq)

    # Case 4: sequence divergence prevents exact match; infer frame from length
    if len(seq) % 3 == 0:
        return "div3_assumed", 0, len(seq)

    return None, None, None


def annotate_alignment(
    gapped_seq: str,
    cds_start: int,
    cds_end: int,
    codon_offset: int = 0,
) -> list[str]:
    """
    Return per-column annotations for the gapped A. thaliana sequence.

    cds_start/cds_end: ungapped positions bounding the coding region.
    codon_offset: for trimmed alignments, the CDS position of the first
                  base (sets the starting frame).

    Annotations: '1'/'2'/'3' = codon position, 'U' = UTR/non-CDS, '-' = gap.
    """
    annotations = []
    ungapped_pos = 0
    codon_counter = codon_offset

    for char in gapped_seq:
        if char in ("-", "."):
            annotations.append("-")
        else:
            if ungapped_pos < cds_start:
                annotations.append("U")
            elif ungapped_pos < cds_end:
                annotations.append(str(codon_counter % 3 + 1))
                codon_counter += 1
            else:
                annotations.append("U")
            ungapped_pos += 1

    return annotations


# ── Diagnostics ───────────────────────────────────────────────────────────────


def run_diagnostics(seq: str, cds: str, cds_start: int, cds_end: int, gene_id: str) -> list[str]:
    """Check CDS sanity. Returns a list of warning strings (empty if all ok)."""
    extracted = seq[cds_start:cds_end]
    starts_atg = extracted[:3] == "ATG"
    stop = extracted[-3:]
    ends_stop = stop in ("TAA", "TAG", "TGA")
    in_frame = len(extracted) % 3 == 0
    cds_length = cds_end - cds_start

    log.info(
        f"{gene_id}: ungapped seq {len(seq)} bp | CDS [{cds_start}:{cds_end}] "
        f"({cds_length} bp) | starts ATG={starts_atg} | ends stop={ends_stop} ({stop}) | "
        f"in-frame={in_frame}"
    )

    warnings = []
    if not starts_atg:
        warnings.append("no_start_ATG")
        log.warning(f"{gene_id}: CDS does not start with ATG")
    if not ends_stop:
        warnings.append(f"no_stop_codon(ends_{stop})")
        log.warning(f"{gene_id}: CDS does not end with stop codon")
    if not in_frame:
        warnings.append(f"not_div3(len={cds_length})")
        log.warning(f"{gene_id}: CDS length not divisible by 3")
    return warnings


# ── Output ────────────────────────────────────────────────────────────────────


def write_annotation(
    out_path: Path,
    annotations: list[str],
    gene_id: str,
    aln_length: int,
):
    """Write a TSV: col_index (0-based), annotation."""
    with open(out_path, "w") as fh:
        fh.write(f"# Gene: {gene_id}\n")
        fh.write(f"# Alignment columns: {aln_length}\n")
        fh.write("# Annotation key: 1/2/3=codon position, U=UTR/non-CDS, -=gap\n")
        fh.write("col_index\tcodon_pos\n")
        for i, ann in enumerate(annotations):
            fh.write(f"{i}\t{ann}\n")
    log.info(f"{gene_id}: annotation written to {out_path}")


# ── Main pipeline ─────────────────────────────────────────────────────────────


def process_gene(
    gene_id: str,
    aln_dir: Path,
    out_dir: Path,
    aln_suffix: str,
    ath_species: str,
    session: requests.Session,
) -> dict:
    """Process a single gene. Returns a status dict for the summary report."""
    result = {
        "gene": gene_id,
        "gene_name": "",
        "status": "ok",
        "confidence": "",
        "ensembl_transcript_id": "",
        "fetched_at_utc": "",
        "note": "",
    }

    # ── Load alignment ────────────────────────────────────────────────────────
    aln_file = aln_dir / f"{gene_id}{aln_suffix}"
    if not aln_file.exists():
        result["status"] = "skip"
        result["note"] = f"alignment file not found: {aln_file}"
        log.warning(result["note"])
        return result

    try:
        alignment = AlignIO.read(str(aln_file), "fasta")
    except Exception as e:
        result["status"] = "error"
        result["note"] = f"could not read alignment: {e}"
        log.error(f"{gene_id}: {result['note']}")
        return result

    aln_length = alignment.get_alignment_length()

    # ── Extract gene name from alignment headers ───────────────────────────────
    try:
        gene_name = extract_gene_name(alignment)
    except ValueError as e:
        result["status"] = "error"
        result["note"] = str(e)
        log.error(f"{gene_id}: {result['note']}")
        return result

    result["gene_name"] = gene_name
    log.info(f"{gene_id}: gene symbol = {gene_name}")

    # ── Find A. thaliana sequence ─────────────────────────────────────────────
    try:
        ath_rec = find_ath_record(alignment, ath_species)
    except ValueError as e:
        result["status"] = "error"
        result["note"] = str(e)
        log.error(f"{gene_id}: {result['note']}")
        return result

    gapped_seq = str(ath_rec.seq).upper()
    ungapped_seq = ungap(gapped_seq)

    # ── Fetch CDS from Ensembl by gene symbol ─────────────────────────────────
    cds, transcript_id = fetch_cds_ensembl(gene_name, session)
    if cds is None:
        result["status"] = "error"
        result["note"] = f"CDS fetch failed for symbol '{gene_name}'"
        return result
    result["ensembl_transcript_id"] = transcript_id
    result["fetched_at_utc"] = datetime.now(UTC).isoformat(timespec="seconds")

    # ── Locate CDS in ungapped sequence ──────────────────────────────────────
    mode, start, end = find_cds_in_seq(ungapped_seq, cds)

    if mode is None:
        result["status"] = "error"
        result["note"] = (
            f"CDS (len={len(cds)}) not found in ungapped A. thaliana seq "
            f"(len={len(ungapped_seq)}); sequences may differ between TAIR versions"
        )
        log.error(f"{gene_id}: {result['note']}")
        return result

    if mode == "flanked":
        cds_start, cds_end, codon_offset = start, end, 0
        result["confidence"] = "high"
        log.info(f"{gene_id}: exact match — CDS found within alignment (has flanking)")
    elif mode == "trimmed":
        cds_start, cds_end, codon_offset = 0, len(ungapped_seq), start
        result["confidence"] = "high"
        log.info(f"{gene_id}: exact match — alignment is CDS subset at offset {start} (frame {start % 3})")
        result["note"] = f"assumed:frame_known trimmed CDS offset={start}"
    elif mode == "near_full":
        cds_start, cds_end, codon_offset = 0, len(ungapped_seq), 0
        diff = len(ungapped_seq) - len(cds)
        result["status"] = "low_confidence"
        result["confidence"] = "low"
        log.info(f"{gene_id}: assumed frame 0 — length differs from CDS by {diff} bp (version mismatch?)")
        result["note"] = f"assumed:frame0 length_diff={diff}bp (version mismatch)"
    else:  # div3_assumed
        cds_start, cds_end, codon_offset = 0, len(ungapped_seq), 0
        result["status"] = "low_confidence"
        result["confidence"] = "low"
        log.info(f"{gene_id}: assumed frame 0 — seq divergence, alignment len {len(ungapped_seq)} divisible by 3")
        result["note"] = f"assumed:frame0 seq_divergence aln={len(ungapped_seq)}bp cds={len(cds)}bp"

    # ── Diagnostics ───────────────────────────────────────────────────────────
    diag_warnings = run_diagnostics(ungapped_seq, cds, cds_start, cds_end, gene_name)
    if diag_warnings:
        result["note"] = (result["note"] + " " if result["note"] else "") + " ".join(diag_warnings)

    # ── Annotate ──────────────────────────────────────────────────────────────
    annotations = annotate_alignment(gapped_seq, cds_start, cds_end, codon_offset)
    assert len(annotations) == aln_length, (
        f"Annotation length {len(annotations)} != alignment length {aln_length}"
    )

    # ── Write output ──────────────────────────────────────────────────────────
    out_path = out_dir / f"{gene_id}_{gene_name}_codon_pos.tsv"
    write_annotation(out_path, annotations, f"{gene_id} ({gene_name})", aln_length)

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", required=True, help="File with one numeric gene ID per line")
    parser.add_argument("--aln_dir", required=True, help="Directory containing alignment files")
    parser.add_argument("--out_dir", required=True, help="Directory for output TSV files")
    parser.add_argument("--aln_suffix", default=".dna.aln.fasta", help="Alignment file suffix (default: .dna.aln.fasta)")
    parser.add_argument(
        "--ath_species",
        default="Arabidopsis_thaliana",
        help="Species string to match in alignment headers (default: Arabidopsis_thaliana)",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.3,
        help="Seconds to wait between Ensembl API calls (default: 0.3)",
    )
    args = parser.parse_args()

    # ── Setup ─────────────────────────────────────────────────────────────────
    aln_dir = Path(args.aln_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_ids = [
        line.strip()
        for line in Path(args.genes).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    log.info(f"Loaded {len(gene_ids)} gene IDs")

    session = requests.Session()
    session.headers.update(HEADERS)

    # ── Process ───────────────────────────────────────────────────────────────
    results = []
    for i, gene_id in enumerate(gene_ids):
        log.info(f"── [{i+1}/{len(gene_ids)}] {gene_id} ──")
        result = process_gene(
            gene_id=gene_id,
            aln_dir=aln_dir,
            out_dir=out_dir,
            aln_suffix=args.aln_suffix,
            ath_species=args.ath_species,
            session=session,
        )
        results.append(result)
        time.sleep(args.delay)

    # ── Summary ───────────────────────────────────────────────────────────────
    summary_path = out_dir / "summary.tsv"
    with open(summary_path, "w") as fh:
        fh.write("gene\tgene_name\tstatus\tconfidence\tensembl_transcript_id\tfetched_at_utc\tnote\n")
        for r in results:
            fh.write(
                f"{r['gene']}\t{r['gene_name']}\t{r['status']}\t{r['confidence']}\t"
                f"{r['ensembl_transcript_id']}\t{r['fetched_at_utc']}\t{r['note']}\n"
            )

    ok = sum(1 for r in results if r["status"] == "ok")
    low_conf = sum(1 for r in results if r["status"] == "low_confidence")
    skip = sum(1 for r in results if r["status"] == "skip")
    err = sum(1 for r in results if r["status"] == "error")
    log.info(f"Done. {ok} ok | {low_conf} low confidence | {skip} skipped | {err} errors — see {summary_path}")


if __name__ == "__main__":
    main()
