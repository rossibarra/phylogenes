#!/usr/bin/env python3
"""
Build derived result files from annotation outputs.

This script classifies summary rows, copies high-confidence FASTA/TSV pairs,
and optionally subsets high-confidence alignments to a species shortlist.
"""

import argparse
import csv
import re
import shutil
from pathlib import Path

from Bio import SeqIO


def read_summary(summary_path: Path) -> list[dict]:
    with open(summary_path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def category_for_row(row: dict) -> str:
    status = row["status"]
    note = row.get("note", "")
    gene_name = row.get("gene_name", "")

    if status == "ok" and note.startswith("assumed:frame_known trimmed"):
        return "exact_trimmed"
    if status == "ok" and (note == "" or row.get("confidence", "") == "high"):
        return "exact_flanked"
    if status in {"ok", "low_confidence"} and note.startswith("assumed:frame0 seq_divergence"):
        return "assumed_seq_divergence"
    if status in {"ok", "low_confidence"} and note.startswith("assumed:frame0 length_diff"):
        return "assumed_version_mismatch"
    if status == "error" and gene_name.startswith("LOC_Os"):
        return "no_arabidopsis_entry"
    if status == "error" and "Cannot find sequence" in note:
        return "no_arabidopsis_sequence"
    if status == "error" and "CDS fetch failed" in note:
        return "cds_fetch_failed"
    if status == "error":
        return "error_other"
    return "other"


def is_high_confidence(row: dict) -> bool:
    note = row.get("note", "")
    confidence = row.get("confidence", "")
    return row["status"] == "ok" and (
        confidence == "high" or note.startswith("assumed:frame_known trimmed") or note == ""
    )


def write_per_gene_notes(rows: list[dict], out_path: Path) -> None:
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_id", "gene_name", "status", "category", "note"])
        for row in rows:
            writer.writerow(
                [
                    row["gene"],
                    row.get("gene_name", ""),
                    row["status"],
                    category_for_row(row),
                    row.get("note", ""),
                ]
            )


def find_annotation_file(annotations_dir: Path, row: dict) -> Path:
    gene = row["gene"]
    gene_name = row.get("gene_name", "")
    expected = annotations_dir / f"{gene}_{gene_name}_codon_pos.tsv"
    if expected.exists():
        return expected
    matches = sorted(annotations_dir.glob(f"{gene}_*_codon_pos.tsv"))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise FileNotFoundError(f"no annotation file found for {gene}")
    raise FileExistsError(f"multiple annotation files found for {gene}: {matches}")


def copy_high_confidence(rows: list[dict], aln_dir: Path, annotations_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        if not is_high_confidence(row):
            continue
        gene = row["gene"]
        aln_file = aln_dir / f"{gene}.dna.aln.fasta"
        annotation_file = find_annotation_file(annotations_dir, row)
        shutil.copy2(aln_file, out_dir / aln_file.name)
        shutil.copy2(annotation_file, out_dir / annotation_file.name)


def species_from_description(description: str) -> str:
    match = re.search(r"Species:(\S+)", description)
    return match.group(1) if match else ""


def normalize_species_name(name: str) -> str:
    return name.strip().replace(" ", "_")


def read_shortlist(shortlist_path: Path, aliases: list[str]) -> list[str]:
    species = [
        normalize_species_name(line)
        for line in shortlist_path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    replacements = {}
    for alias in aliases:
        source, target = alias.split("=", 1)
        replacements[normalize_species_name(source)] = normalize_species_name(target)
    return [replacements.get(item, item) for item in species]


def record_matches_species(record_species: str, shortlist_species: list[str]) -> bool:
    return any(record_species == species or record_species.startswith(f"{species}_") for species in shortlist_species)


def subset_shortlist_alignments(
    rows: list[dict],
    high_confidence_dir: Path,
    shortlist_path: Path,
    out_dir: Path,
    aliases: list[str],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    shortlist_species = read_shortlist(shortlist_path, aliases)
    for row in rows:
        if not is_high_confidence(row):
            continue
        gene = row["gene"]
        source_fasta = high_confidence_dir / f"{gene}.dna.aln.fasta"
        records = [
            rec for rec in SeqIO.parse(source_fasta, "fasta")
            if record_matches_species(species_from_description(rec.description), shortlist_species)
        ]
        SeqIO.write(records, out_dir / source_fasta.name, "fasta")
        annotation_file = find_annotation_file(high_confidence_dir, row)
        shutil.copy2(annotation_file, out_dir / annotation_file.name)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", default="results/codon_annotations/summary.tsv")
    parser.add_argument("--aln_dir", default="data/alignments")
    parser.add_argument("--annotations_dir", default="results/codon_annotations")
    parser.add_argument("--per_gene_notes", default="results/per_gene_notes.tsv")
    parser.add_argument("--high_confidence_dir", default="results/high_confidence_alignments")
    parser.add_argument("--shortlist", default="shortlist.txt")
    parser.add_argument("--shortlist_dir", default="results/shortlist_files")
    parser.add_argument(
        "--species_alias",
        action="append",
        default=[],
        help="Map one shortlist species to another, e.g. 'Streptanthus tortuosus=Streptanthus carinatus'",
    )
    parser.add_argument("--skip_shortlist", action="store_true")
    args = parser.parse_args()

    rows = read_summary(Path(args.summary))
    write_per_gene_notes(rows, Path(args.per_gene_notes))
    copy_high_confidence(
        rows,
        Path(args.aln_dir),
        Path(args.annotations_dir),
        Path(args.high_confidence_dir),
    )
    if not args.skip_shortlist:
        subset_shortlist_alignments(
            rows,
            Path(args.high_confidence_dir),
            Path(args.shortlist),
            Path(args.shortlist_dir),
            args.species_alias,
        )


if __name__ == "__main__":
    main()
