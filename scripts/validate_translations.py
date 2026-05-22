#!/usr/bin/env python3
"""
validate_translations.py

For a random sample of high-confidence annotated genes:
  1. Extracts the A. thaliana CDS from the alignment using the codon annotation
  2. Translates it
  3. Fetches the canonical protein from Ensembl
  4. Reports identity and any mismatches
"""

import argparse
import random
import re
import sys
from pathlib import Path

import requests
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq

ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}


def fetch_protein_ensembl(gene_symbol: str, session: requests.Session) -> tuple[str, str] | tuple[None, None]:
    """Return (transcript_id, protein_seq) for the canonical A. thaliana transcript."""
    if re.match(r'^AT[1-5MC]G\d+$', gene_symbol, re.IGNORECASE):
        url = f"{ENSEMBL_REST}/lookup/id/{gene_symbol}"
    else:
        url = f"{ENSEMBL_REST}/lookup/symbol/arabidopsis_thaliana/{gene_symbol}"

    try:
        r = session.get(url, headers=HEADERS, params={"expand": 1}, timeout=30)
        r.raise_for_status()
        gene_data = r.json()
    except Exception as e:
        return None, f"gene lookup failed: {e}"

    transcripts = gene_data.get("Transcript", [])
    if not transcripts:
        return None, "no transcripts"

    canonical = next((t for t in transcripts if t.get("is_canonical") == 1), None)
    if canonical is None:
        canonical = max(
            transcripts,
            key=lambda t: t.get("Translation", {}).get("length", 0) if t.get("Translation") else 0,
        )

    transcript_id = canonical["id"]
    if not canonical.get("Translation"):
        return transcript_id, "no translation (non-coding?)"

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
        return transcript_id, f"protein fetch failed: {e}"

    if isinstance(data, list):
        data = data[0]
    protein = data.get("seq", "").upper()
    return transcript_id, protein


def extract_cds_from_alignment(aln_file: Path, codon_offset: int, ath_species: str = "Arabidopsis_thaliana") -> str | None:
    """Extract the ungapped A. thaliana sequence from the alignment."""
    alignment = AlignIO.read(str(aln_file), "fasta")
    ath_rec = next((r for r in alignment if ath_species in r.description), None)
    if ath_rec is None:
        return None
    return str(ath_rec.seq).replace("-", "").replace(".", "").upper()


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
    expected_aa_start = codon_offset // 3

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
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    # Load high-confidence trimmed genes with their offsets
    genes = []
    with open(args.summary) as f:
        next(f)
        for line in f:
            gene, gene_name, status, note = line.rstrip("\n").split("\t")
            if status == "ok" and note.startswith("assumed:frame_known trimmed"):
                m = re.search(r"offset=(\d+)", note)
                offset = int(m.group(1)) if m else 0
                genes.append((gene, gene_name, offset))

    random.seed(args.seed)
    sample = random.sample(genes, min(args.n, len(genes)))

    session = requests.Session()
    session.headers.update(HEADERS)

    hc_dir = Path(args.hc_dir)

    print(f"{'gene':<8} {'symbol':<16} {'offset':>6}  {'match':<8} {'identity':>8}  detail")
    print("-" * 90)

    for gene, gene_name, offset in sample:
        aln_file = hc_dir / f"{gene}.dna.aln.fasta"
        seq = extract_cds_from_alignment(aln_file, offset)
        if seq is None:
            print(f"{gene:<8} {gene_name:<16} {offset:>6}  ERROR    —  no A. thaliana sequence")
            continue

        our_prot = translate_cds(seq, offset)
        transcript_id, ref_prot = fetch_protein_ensembl(gene_name, session)

        if ref_prot is None or not ref_prot.replace("*", "").isalpha():
            print(f"{gene:<8} {gene_name:<16} {offset:>6}  ERROR    —  {ref_prot}")
            continue

        result = compare_proteins(our_prot, ref_prot, offset)
        print(
            f"{gene:<8} {gene_name:<16} {offset:>6}  {result['match']:<8} {result['identity']:>8.1%}  {result['detail']}"
        )

        if result["match"] not in ("exact", "good"):
            aa_start = offset // 3
            ref_slice = ref_prot.rstrip("*")[aa_start: aa_start + len(our_prot)]
            print(f"         ours: {our_prot[:80]}{'...' if len(our_prot)>80 else ''}")
            print(f"         ref:  {ref_slice[:80]}{'...' if len(ref_slice)>80 else ''}")


if __name__ == "__main__":
    main()
