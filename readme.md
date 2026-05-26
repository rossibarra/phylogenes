# phylogenes

Annotations of multi-species gene alignments in the context of *Arabidopsis thaliana* gene models.

## Overview

This repo annotates codon positions across large phylogenomic DNA alignments sourced from the [Kew Tree of Life Explorer](https://sftp.kew.org/pub/treeoflife/current_release/fasta/alignments/). Each alignment covers hundreds of species for a single gene family. The pipeline fetches the canonical *A. thaliana* CDS from Ensembl Plants, locates the *A. thaliana* sequence within each alignment, and annotates every alignment column with its codon position (1, 2, or 3) or gap status.

## Input data

Alignment files are named `{numeric_id}.dna.aln.fasta` and downloaded from the Kew Tree of Life Explorer:

```
https://sftp.kew.org/pub/treeoflife/current_release/fasta/alignments/
```

353 alignments are stored in `data/alignments/` (listed with sizes in `data/alignments/file_list.txt`). Files are excluded from version control (see `.gitignore`) due to size (~6 GB total). To re-download:

```bash
awk '{print "https://sftp.kew.org/pub/treeoflife/current_release/fasta/alignments/" $1}' \
    data/alignments/file_list.txt | wget -i - -P data/alignments/ -c
```

FASTA headers follow this format:

```
>{id} Gene_Name:{symbol} Species:{species} Repository:INSDC Sequence_ID:{accession}
```

`data/gene_list.txt` contains one numeric ID per line. `data/alignments/file_list.txt` lists the available alignment files with sizes.

## Scripts

### `scripts/annotate_codon_positions.py`

Main annotation pipeline. For each gene:

1. Loads the alignment (`{id}.dna.aln.fasta`)
2. Extracts the gene symbol from the `Gene_Name:` field in the alignment headers
3. Fetches the canonical CDS from the [Ensembl Plants REST API](https://rest.ensembl.org) using the gene symbol (or TAIR ID directly for `AT[1-5]G`-format names)
4. Finds the *A. thaliana* sequence by matching `Species:Arabidopsis_thaliana` in the header description
5. Locates the alignment sequence within the CDS to determine the reading frame offset
6. Annotates every alignment column with its codon position
7. Writes a TSV annotation file per gene and a `summary.tsv`

**Usage:**
```bash
python scripts/annotate_codon_positions.py \
    --genes data/gene_list.txt \
    --aln_dir data/alignments \
    --out_dir results/codon_annotations \
    [--aln_suffix .dna.aln.fasta] \
    [--ath_species Arabidopsis_thaliana] \
    [--delay 0.3]
```

**Output:** one TSV per gene (`{id}_{symbol}_codon_pos.tsv`) plus `summary.tsv`.

| Column | Values |
|--------|--------|
| `col_index` | 0-based alignment column index |
| `codon_pos` | `1`, `2`, `3` = codon position; `U` = UTR/non-CDS; `-` = gap column |

**`summary.tsv` columns:** `gene`, `gene_name`, `status`, `confidence`, `ensembl_transcript_id`, `fetched_at_utc`, `note`

`status=ok` is reserved for exact CDS matches with high-confidence frame annotations. Records where the exact CDS position could not be established are written as `status=low_confidence` with `confidence=low`.

The `note` field encodes how the frame was determined and any diagnostic warnings (space-separated):

| Note prefix | Meaning |
|-------------|---------|
| *(empty)* | Exact match: full CDS found within the alignment (alignment has flanking sequence) |
| `assumed:frame_known trimmed CDS offset=N` | Exact match: alignment is a substring of the CDS starting at CDS position N; frame derived from offset |
| `assumed:frame0 length_diff=Nbp (version mismatch)` | Low confidence: lengths differ by ≤9 bp; exact match failed; frame 0 assumed |
| `assumed:frame0 seq_divergence aln=Xbp cds=Ybp` | Low confidence: exact match failed; alignment length divisible by 3; frame 0 assumed |
| `no_start_ATG` | Extracted CDS region does not begin with ATG (expected for trimmed alignments) |
| `no_stop_codon(ends_NNN)` | Extracted CDS region does not end with a stop codon (expected for trimmed alignments) |
| `not_div3(len=N)` | CDS region length not divisible by 3 |

### `scripts/validate_translations.py`

Validation script. For a random sample of annotated genes, reconstructs the annotated *A. thaliana* CDS from the `*_codon_pos.tsv` columns, translates it, and compares the result against the Ensembl reference protein. Checks that the translated slice matches the expected position in the reference.

**Usage:**
```bash
python scripts/validate_translations.py \
    --summary results/codon_annotations/summary.tsv \
    --hc_dir results/high_confidence_alignments \
    --annotation_dir results/high_confidence_alignments \
    --n 10 \
    [--seed 42] \
    [--all]
```

The default `--hc_dir` requires the precomputed high-confidence FASTA files. On a fresh clone, first download `data/alignments/` and run the derived-output command below.

### `scripts/build_derived_outputs.py`

Builds the cleaned category table, high-confidence FASTA/TSV directory, and shortlist-subsetted files from `summary.tsv`.

**Usage:**
```bash
python scripts/build_derived_outputs.py \
    --summary results/codon_annotations/summary.tsv \
    --aln_dir data/alignments \
    --annotations_dir results/codon_annotations \
    --per_gene_notes results/per_gene_notes.tsv \
    --high_confidence_dir results/high_confidence_alignments \
    --shortlist shortlist.txt \
    --shortlist_dir results/shortlist_files \
    --species_alias "Streptanthus tortuosus=Streptanthus carinatus"
```

## Results

### Annotation run (353 alignments, Kew Tree of Life current release)

| Category | Count |
|----------|-------|
| Annotated — exact match, frame known from CDS offset | 222 |
| Annotated — exact/flanked match | 5 |
| Low confidence — frame assumed from sequence divergence | 30 |
| Low confidence — frame assumed from version mismatch | 1 |
| Error — gene symbol not found in Ensembl | 72 |
| Error — rice `LOC_Os` symbol (no *A. thaliana* entry) | 21 |
| Error — no *A. thaliana* sequence in alignment | 1 |
| Error — other | 1 |

Full per-gene results: `results/codon_annotations/summary.tsv`. A cleaned version with standardised category labels is at `results/per_gene_notes.tsv`.

### Validation

The validation script now reconstructs the coding slice from the codon-position TSV before translation, so it checks the annotation file rather than only ungapping the FASTA row. After the confidence metadata update, 10 high-confidence genes were sampled at random (seed 42):

| Gene ID | Symbol | CDS offset | Transcript | Ref slice | Result |
|---------|--------|-----------:|------------|-----------|--------|
| 6631 | AT1G74530 | 177 | AT1G74530.3 | ref[59:309] | exact, 100% |
| 5354 | AT1G18340 | 30 | AT1G18340.1 | ref[10:277] | exact, 100% |
| 4848 | ARC6 | 246 | AT5G42480.1 | ref[82:799] | exact, 100% |
| 6865 | AT1G22800 | 150 | AT1G22800.1 | ref[50:351] | exact, 100% |
| 5843 | AT5G52980 | 207 | AT5G52980.1 | ref[69:187] | exact, 100% |
| 5744 | GDC1 | 171 | AT1G50900.1 | ref[57:173] | exact, 100% |
| 5670 | APC7 | 9 | AT2G39090.1 | ref[3:558] | exact, 100% |
| 5428 | AT1G62730 | 30 | AT1G62730.1 | ref[10:304] | exact, 100% |
| 6864 | NPU | 96 | AT3G51610.1 | ref[32:217] | exact, 100% |
| 5335 | AT4G28450 | 0 | AT4G28450.1 | ref[0:452] | exact, 100% |

**10/10 exact matches.** All TSV-reconstructed translations matched the corresponding slice of the Ensembl reference protein with 100% identity and zero coding gaps. For final analyses, run with `--all` to validate all high-confidence genes.

### High-confidence output

The 227 high-confidence exact-match genes are collected in `results/high_confidence_alignments/`, with both the source alignment and annotation TSV for each gene:

```
results/high_confidence_alignments/
├── 4471.dna.aln.fasta
├── 4471_D2HGDH_codon_pos.tsv
├── 4527.dna.aln.fasta
├── 4527_PUR2_codon_pos.tsv
└── ...  (227 pairs = 454 files total)
```

### Species shortlist and subsetted alignments

`shortlist.txt` contains 16 focal crop and model species. The source list names *Streptanthus tortuosus*; generated shortlist outputs use *Streptanthus carinatus* as the available proxy via the alias shown in the derived-output command above.

| Species | Notes |
|---------|-------|
| *Arabidopsis thaliana* | |
| *Brassica napus* | |
| *Streptanthus tortuosus* | generated outputs use *S. carinatus* as proxy |
| *Phaseolus vulgaris* | |
| *Glycine max* | |
| *Helianthus annuus* | |
| *Capsicum annuum* | |
| *Solanum lycopersicum* | |
| *Amaranthus hypochondriacus* | |
| *Gossypium hirsutum* | |
| *Oryza sativa* | matched on prefix (includes subspecies) |
| *Triticum aestivum* | |
| *Hordeum vulgare* | |
| *Zea mays* | |
| *Setaria italica* | |
| *Sorghum bicolor* | |

For each of the 227 high-confidence genes, alignments were subsetted to these species only and copied alongside their annotation TSVs to `results/shortlist_files/`. Current shortlist FASTAs contain 9 files with 16 records, 36 with 17 records, and 182 with 18 records. Counts above 16 occur because prefix matching includes subspecies such as `Glycine_max_subsp._soja` and `Sorghum_bicolor_nothosubsp._drummondii`. The annotation TSVs are copied as-is since column indices are independent of which species are present.

`results/shortlist_files/gap_stripped/` contains versions of these files with alignment columns that are entirely gaps across all sequences removed. The annotation TSVs are updated accordingly: rows annotated `-` (thaliana gap) that corresponded to all-gap columns are dropped, and `col_index` is renumbered from 0 to match the new column positions. Columns annotated `-` where at least one non-thaliana sequence has a nucleotide are retained. Generated by `scripts/strip_gap_columns.py`.

## Setup

```bash
conda create -n phylogenes "biopython>=1.81" "requests>=2.31"
conda activate phylogenes
```

Or with pip:
```bash
pip install -r requirements.txt
```

## Directory structure

```
phylogenes/
├── data/
│   ├── gene_list.txt              # one numeric gene ID per line
│   └── alignments/
│       ├── file_list.txt          # source file list with sizes
│       └── {id}.dna.aln.fasta    # one alignment per gene
├── shortlist.txt                  # 16 focal species for subsetting
├── results/
│   ├── codon_annotations/         # per-gene TSVs + summary.tsv
│   ├── per_gene_notes.tsv         # cleaned annotation category table
│   ├── high_confidence_alignments/ # exact-match genes: FASTA + TSV pairs
│   └── shortlist_files/           # shortlist-subsetted alignments + TSVs
├── scripts/
│   ├── annotate_codon_positions.py
│   ├── build_derived_outputs.py
│   └── validate_translations.py
├── requirements.txt
└── readme.md
```
