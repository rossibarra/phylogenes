#!/usr/bin/env python3
"""
validate_translations.py

For a random sample of high-confidence annotated genes:
  1. Reconstructs the A. thaliana CDS slice from the alignment and codon annotation
  2. Translates it
  3. Fetches the canonical protein from Ensembl
  4. Reports identity and any mismatches
"""

import argparse
import csv
import random
import re
from pathlib import Path

import requests
from Bio import AlignIO
from Bio.Seq import Seq

ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}


def fetch_protein_ensembl(gene_symbol: str, session: requests.Session) -> tuple[str | None, str | None, str | None]:
    """Return (transcript_id, protein_id, protein_seq) for the canonical A. thaliana transcript."""
    if re.match(r'^AT[1-5MC]G\d+$', gene_symbol, re.IGNORECASE):
        url = f"{ENSEMBL_REST}/lookup/id/{gene_symbol}"
    else:
        url = f"{ENSEMBL_REST}/lookup/symbol/arabidopsis_thaliana/{gene_symbol}"

    try:
        r = session.get(url, headers=HEADERS, params={"expand": 1}, timeout=30)
        r.raise_for_status()
        gene_data = r.json()
    except Exception as e:
        return None, None, f"gene lookup failed: {e}"

    transcripts = gene_data.get("Transcript", [])
    if not transcripts:
        return None, None, "no transcripts"

    canonical = next((t for t in transcripts if t.get("is_canonical") == 1), None)
    if canonical is None:
        canonical = max(
            transcripts,
            key=lambda t: t.get("Translation", {}).get("length", 0) if t.get("Translation") else 0,
        )

    transcript_id = canonical["id"]
    translation = canonical.get("Translation")
    if not translation:
        return transcript_id, None, "no translation (non-coding?)"
    protein_id = translation.get("id")

    try:
        r = session.get(
            f"{ENSEMBL_REST}/sequence/id/{transcript_id}",
            headers=HEADERS,
            params={"type": "protein"},
            timeout=30,
        )
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        return transcript_id, protein_id, f"protein fetch failed: {e}"

    if isinstance(data, list):
        data = data[0]
    protein = data.get("seq", "").upper()
    return transcript_id, protein_id, protein


def read_codon_annotations(annotation_file: Path) -> list[str]:
    """Read the per-column codon-position annotations, ignoring metadata comments."""
    annotations = []
    with open(annotation_file, newline="") as fh:
        rows = (line for line in fh if not line.startswith("#"))
        reader = csv.DictReader(rows, delimiter="\t")
        if reader.fieldnames != ["col_index", "codon_pos"]:
            raise ValueError(f"{annotation_file}: expected columns col_index,codon_pos")
        for expected_index, row in enumerate(reader):
            col_index = int(row["col_index"])
            if col_index != expected_index:
                raise ValueError(
                    f"{annotation_file}: non-contiguous col_index at row {expected_index}: {col_index}"
                )
            annotations.append(row["codon_pos"])
    return annotations


def find_annotation_file(annotation_dir: Path, gene: str, gene_name: str) -> Path:
    """Find the annotation TSV for a gene, tolerating changed gene-symbol casing."""
    expected = annotation_dir / f"{gene}_{gene_name}_codon_pos.tsv"
    if expected.exists():
        return expected
    matches = sorted(annotation_dir.glob(f"{gene}_*_codon_pos.tsv"))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise FileNotFoundError(f"no annotation file found for gene {gene} in {annotation_dir}")
    raise FileExistsError(f"multiple annotation files found for gene {gene}: {matches}")


def extract_annotated_cds_from_alignment(
    aln_file: Path,
    annotation_file: Path,
    ath_species: str = "Arabidopsis_thaliana",
) -> tuple[str | None, dict]:
    """Reconstruct the annotated A. thaliana coding slice from FASTA columns labeled 1/2/3."""
    alignment = AlignIO.read(str(aln_file), "fasta")
    ath_rec = next((r for r in alignment if ath_species in r.description), None)
    if ath_rec is None:
        return None, {"error": "no A. thaliana sequence"}

    gapped_seq = str(ath_rec.seq).upper()
    annotations = read_codon_annotations(annotation_file)
    if len(annotations) != len(gapped_seq):
        raise ValueError(
            f"{annotation_file}: {len(annotations)} annotations for {len(gapped_seq)} alignment columns"
        )

    bases = []
    coding_positions = 0
    coding_gaps = 0
    phase_errors = []
    for col_index, (base, codon_pos) in enumerate(zip(gapped_seq, annotations, strict=True)):
        if codon_pos in {"1", "2", "3"}:
            coding_positions += 1
            expected = str((len(bases) % 3) + 1)
            if codon_pos != expected:
                phase_errors.append((col_index, codon_pos, expected))
            if base in {"-", "."}:
                coding_gaps += 1
            else:
                bases.append(base)
        elif codon_pos == "-":
            if base not in {"-", "."}:
                phase_errors.append((col_index, codon_pos, "base"))
        elif codon_pos != "U":
            raise ValueError(f"{annotation_file}: invalid codon_pos '{codon_pos}' at column {col_index}")

    diagnostics = {
        "coding_columns": coding_positions,
        "coding_gaps": coding_gaps,
        "phase_errors": phase_errors,
        "annotation_file": str(annotation_file),
    }
    return "".join(bases), diagnostics


def translate_cds(seq: str, codon_offset: int) -> str:
    """
    Translate seq given that it starts at CDS position `codon_offset`.
    Skips incomplete leading codon, translates all complete in-frame codons.
    Strips trailing stop if present.
    """
    skip = (3 - codon_offset % 3) % 3
    coding = seq[skip:]
    # trim to multiple of 3
    coding = coding[: len(coding) - len(coding) % 3]
    protein = str(Seq(coding).translate())
    # strip trailing stop
    if protein.endswith("*"):
        protein = protein[:-1]
    return protein


def translated_ref_start(codon_offset: int) -> int:
    """Reference amino-acid start after dropping an incomplete leading codon."""
    skip = (3 - codon_offset % 3) % 3
    return (codon_offset + skip) // 3


def compare_proteins(our_prot: str, ref_prot: str, codon_offset: int) -> dict:
    """
    Compare our translated slice against the reference protein.
    Since our translation covers a trimmed region of the CDS, we check
    whether our protein matches the expected slice of the reference:
      ref[offset//3 : offset//3 + len(ours)]
    Also checks if our protein is present anywhere as a substring (frame sanity).
    """
    ref = ref_prot.rstrip("*")
    our = our_prot
    expected_aa_start = translated_ref_start(codon_offset)

    # Primary check: does our protein match the expected position in the reference?
    ref_slice = ref[expected_aa_start: expected_aa_start + len(our)]
    if ref_slice == our:
        return {
            "match": "exact",
            "identity": 1.0,
            "detail": f"matches ref[{expected_aa_start}:{expected_aa_start+len(our)}] perfectly",
        }

    if len(ref_slice) == 0:
        return {"match": "out_of_range", "identity": 0.0,
                "detail": f"expected_aa_start={expected_aa_start} beyond ref len={len(ref)}"}

    # Identity at expected position
    matches = sum(a == b for a, b in zip(our, ref_slice))
    identity = matches / max(len(our), len(ref_slice))
    mismatches = [(i, our[i], ref_slice[i]) for i in range(min(len(our), len(ref_slice))) if our[i] != ref_slice[i]]
    detail = (f"ref[{expected_aa_start}:{expected_aa_start+len(our)}] "
              f"len_ours={len(our)} len_ref_slice={len(ref_slice)} "
              f"mismatches={len(mismatches)}")
    if mismatches[:3]:
        detail += " eg:" + ",".join(f"{i}:{o}>{r}" for i, o, r in mismatches[:3])

    # Substring check: does our protein appear anywhere in the reference?
    sub_pos = ref.find(our[:20]) if len(our) >= 20 else -1
    if sub_pos != -1:
        detail += f" (found in ref at aa {sub_pos})"

    return {
        "match": "good" if identity >= 0.95 else "partial" if identity >= 0.5 else "poor",
        "identity": identity,
        "detail": detail,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", default="results/codon_annotations/summary.tsv")
    parser.add_argument("--hc_dir", default="results/high_confidence_alignments")
    parser.add_argument(
        "--annotation_dir",
        default=None,
        help="Directory containing *_codon_pos.tsv files (default: --hc_dir)",
    )
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--all", action="store_true", help="Validate all high-confidence genes")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    # Load high-confidence trimmed genes with their offsets
    genes = []
    with open(args.summary) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row["gene"]
            gene_name = row["gene_name"]
            status = row["status"]
            confidence = row.get("confidence", "")
            note = row["note"]
            is_high_confidence = (
                status == "ok"
                and (
                    confidence == "high"
                    or note.startswith("assumed:frame_known trimmed")
                    or note == ""
                )
            )
            if is_high_confidence and note.startswith("assumed:frame_known trimmed"):
                m = re.search(r"offset=(\d+)", note)
                offset = int(m.group(1)) if m else 0
                genes.append((gene, gene_name, offset))
            elif is_high_confidence and confidence == "high":
                genes.append((gene, gene_name, 0))

    if args.all:
        sample = genes
    else:
        random.seed(args.seed)
        sample = random.sample(genes, min(args.n, len(genes)))

    session = requests.Session()
    session.headers.update(HEADERS)

    hc_dir = Path(args.hc_dir)
    annotation_dir = Path(args.annotation_dir) if args.annotation_dir else hc_dir

    print(
        f"{'gene':<8} {'symbol':<16} {'offset':>6}  {'match':<8} {'identity':>8}  "
        f"{'transcript':<18} {'protein':<18} detail"
    )
    print("-" * 130)

    for gene, gene_name, offset in sample:
        aln_file = hc_dir / f"{gene}.dna.aln.fasta"
        try:
            annotation_file = find_annotation_file(annotation_dir, gene, gene_name)
            seq, diagnostics = extract_annotated_cds_from_alignment(aln_file, annotation_file)
        except Exception as e:
            print(f"{gene:<8} {gene_name:<16} {offset:>6}  ERROR    —  {'':<18} {'':<18} {e}")
            continue

        if seq is None:
            print(f"{gene:<8} {gene_name:<16} {offset:>6}  ERROR    —  {'':<18} {'':<18} {diagnostics['error']}")
            continue
        if diagnostics["phase_errors"]:
            print(
                f"{gene:<8} {gene_name:<16} {offset:>6}  ERROR    —  {'':<18} {'':<18} "
                f"annotation phase errors: {diagnostics['phase_errors'][:3]}"
            )
            continue

        our_prot = translate_cds(seq, offset)
        transcript_id, protein_id, ref_prot = fetch_protein_ensembl(gene_name, session)

        if ref_prot is None or not ref_prot.replace("*", "").isalpha():
            print(
                f"{gene:<8} {gene_name:<16} {offset:>6}  ERROR    —  "
                f"{transcript_id or '':<18} {protein_id or '':<18} {ref_prot}"
            )
            continue

        result = compare_proteins(our_prot, ref_prot, offset)
        detail = (
            f"{result['detail']} | tsv_cols={diagnostics['coding_columns']} "
            f"coding_gaps={diagnostics['coding_gaps']}"
        )
        print(
            f"{gene:<8} {gene_name:<16} {offset:>6}  {result['match']:<8} {result['identity']:>8.1%}  "
            f"{transcript_id or '':<18} {protein_id or '':<18} {detail}"
        )

        if result["match"] not in ("exact", "good"):
            aa_start = offset // 3
            ref_slice = ref_prot.rstrip("*")[aa_start: aa_start + len(our_prot)]
            print(f"         ours: {our_prot[:80]}{'...' if len(our_prot)>80 else ''}")
            print(f"         ref:  {ref_slice[:80]}{'...' if len(ref_slice)>80 else ''}")


if __name__ == "__main__":
    main()
