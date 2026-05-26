#!/usr/bin/env python3
"""
strip_gap_columns.py

For each gene alignment in a source directory:
  1. Removes alignment columns that are entirely gaps across all sequences.
  2. Writes the stripped alignment to the output directory.
  3. Rewrites the codon_pos TSV: drops '-' rows and renumbers col_index.

Usage:
    python scripts/strip_gap_columns.py \
        --in_dir  results/shortlist_files \
        --out_dir results/gap_stripped \
        [--aln_suffix .dna.aln.fasta]
"""

import argparse
import logging
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


def gap_only_columns(alignment) -> set[int]:
    """Return the set of 0-based column indices that are all gaps."""
    n_cols = alignment.get_alignment_length()
    gap_cols = set()
    for col in range(n_cols):
        if all(rec.seq[col] == "-" for rec in alignment):
            gap_cols.add(col)
    return gap_cols


def strip_columns(alignment, drop: set[int]) -> MultipleSeqAlignment:
    """Return a new alignment with the specified columns removed."""
    n_cols = alignment.get_alignment_length()
    keep = [i for i in range(n_cols) if i not in drop]
    new_records = []
    for rec in alignment:
        new_seq = "".join(rec.seq[i] for i in keep)
        new_records.append(
            SeqRecord(Seq(new_seq), id=rec.id, description=rec.description)
        )
    return MultipleSeqAlignment(new_records)


def rewrite_tsv(in_tsv: Path, out_tsv: Path, drop: set[int]) -> None:
    """
    Copy the codon_pos TSV, skipping rows whose col_index is in drop,
    and renumbering the remaining col_index values sequentially from 0.
    Comment lines (starting with #) are updated where they mention column counts.
    """
    lines = in_tsv.read_text().splitlines()
    header_comments = []
    data_rows = []  # (original_col_index, codon_pos)

    for line in lines:
        if line.startswith("#"):
            header_comments.append(line)
        elif line.startswith("col_index"):
            pass  # skip column header; we'll re-emit it
        elif line.strip():
            parts = line.split("\t")
            orig_idx = int(parts[0])
            cpos = parts[1]
            data_rows.append((orig_idx, cpos))

    # Patch the "Alignment columns: N" comment to reflect new count
    kept_count = sum(1 for idx, _ in data_rows if idx not in drop)
    new_comments = []
    for c in header_comments:
        if "Alignment columns:" in c:
            c = c.split("Alignment columns:")[0] + f"Alignment columns: {kept_count}"
        new_comments.append(c)

    out_lines = new_comments + ["col_index\tcodon_pos"]
    new_idx = 0
    for orig_idx, cpos in data_rows:
        if orig_idx in drop:
            continue
        out_lines.append(f"{new_idx}\t{cpos}")
        new_idx += 1

    out_tsv.write_text("\n".join(out_lines) + "\n")


def process_gene(aln_path: Path, tsv_path: Path, out_dir: Path, aln_suffix: str) -> None:
    alignment = AlignIO.read(str(aln_path), "fasta")
    drop = gap_only_columns(alignment)

    stripped = strip_columns(alignment, drop)
    out_aln = out_dir / aln_path.name
    with open(out_aln, "w") as fh:
        for rec in stripped:
            fh.write(f">{rec.description}\n{rec.seq}\n")

    out_tsv = out_dir / tsv_path.name
    rewrite_tsv(tsv_path, out_tsv, drop)

    log.info(
        "%s: %d columns → %d kept, %d gap-only dropped",
        aln_path.name,
        alignment.get_alignment_length(),
        alignment.get_alignment_length() - len(drop),
        len(drop),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in_dir",  required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--aln_suffix", default=".dna.aln.fasta")
    args = parser.parse_args()

    in_dir  = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    aln_files = sorted(in_dir.glob(f"*{args.aln_suffix}"))
    if not aln_files:
        log.error("No alignment files found in %s", in_dir)
        return

    for aln_path in aln_files:
        # Find the matching TSV: same numeric prefix, ends with _codon_pos.tsv
        prefix = aln_path.name.replace(args.aln_suffix, "")
        tsv_candidates = sorted(in_dir.glob(f"{prefix}_*_codon_pos.tsv"))
        if not tsv_candidates:
            log.warning("No TSV found for %s — skipping", aln_path.name)
            continue
        process_gene(aln_path, tsv_candidates[0], out_dir, args.aln_suffix)


if __name__ == "__main__":
    main()
