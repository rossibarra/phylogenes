# phylogenes

Tools for annotating and analyzing multi-species gene alignments in the context of *Arabidopsis thaliana* gene models.

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

**`summary.tsv` columns:** `gene`, `gene_name`, `status`, `note`

The `note` field encodes how the frame was determined and any diagnostic warnings (space-separated):

| Note prefix | Meaning |
|-------------|---------|
| *(empty)* | Exact match: full CDS found within the alignment (alignment has flanking sequence) |
| `assumed:frame_known trimmed CDS offset=N` | Exact match: alignment is a substring of the CDS starting at CDS position N; frame derived from offset |
| `assumed:frame0 length_diff=Nbp (version mismatch)` | Lengths differ by ≤9 bp; exact match failed; frame 0 assumed |
| `assumed:frame0 seq_divergence aln=Xbp cds=Ybp` | Exact match failed; alignment length divisible by 3; frame 0 assumed |
| `no_start_ATG` | Extracted CDS region does not begin with ATG (expected for trimmed alignments) |
| `no_stop_codon(ends_NNN)` | Extracted CDS region does not end with a stop codon (expected for trimmed alignments) |
| `not_div3(len=N)` | CDS region length not divisible by 3 |

### `scripts/validate_translations.py`

Validation script. For a random sample of annotated genes, translates the annotated *A. thaliana* CDS and compares the result against the Ensembl reference protein. Checks that the translated slice matches the expected position in the reference (`ref[offset//3 : offset//3 + len(translation)]`).

**Usage:**
```bash
python scripts/validate_translations.py \
    --summary results/codon_annotations/summary.tsv \
    --hc_dir results/high_confidence_alignments \
    --n 10 \
    [--seed 42]
```

## Results

### Annotation run (353 alignments, Kew Tree of Life current release)

| Category | Count |
|----------|-------|
| Annotated — exact match, frame known from CDS offset | 222 |
| Annotated — frame assumed (seq divergence or version mismatch) | 31 |
| Annotated — exact match with flanking sequence | ~5 |
| Error — gene symbol not found in Ensembl | 72 |
| Error — rice `LOC_Os` symbol (no *A. thaliana* entry) | 18 |
| Error — no *A. thaliana* sequence in alignment | 4 |
| Error — other | 1 |

Full per-gene results: `results/codon_annotations/summary.tsv`. A cleaned version with standardised category labels is at `results/per_gene_notes.tsv`.

### Validation

The 222 high-confidence annotations were validated by translating each annotated *A. thaliana* region and comparing the result against the Ensembl reference protein at the expected position (`ref[offset//3 : offset//3 + len(translation)]`). 10 genes were sampled at random (seed 42):

| Gene ID | Symbol | CDS offset | Ref slice | Result |
|---------|--------|-----------|-----------|--------|
| 6652 | EBS3 | 144 | ref[48:567] | ✓ exact |
| 5355 | PyrD | 216 | ref[72:387] | ✓ exact |
| 4848 | ARC6 | 246 | ref[82:799] | ✓ exact |
| 6914 | AT5G06260 | 33 | ref[11:422] | ✓ exact |
| 5858 | AT2G40650 | 0 | ref[0:188] | ✓ exact |
| 5772 | PDF1A | 240 | ref[80:267] | ✓ exact |
| 5703 | EMB2762 | 60 | ref[20:568] | ✓ exact |
| 5430 | AT4G29520 | 102 | ref[34:217] | ✓ exact |
| 6913 | HCF136 | 330 | ref[110:426] | ✓ exact |
| 5348 | SLD5 | 54 | ref[18:220] | ✓ exact |

**10/10 exact matches.** All translated sequences matched the corresponding slice of the Ensembl reference protein with 100% identity.

### High-confidence output

The 222 exact-match genes are collected in `results/high_confidence_alignments/`, with both the source alignment and annotation TSV for each gene:

```
results/high_confidence_alignments/
├── 4471.dna.aln.fasta
├── 4471_D2HGDH_codon_pos.tsv
├── 4527.dna.aln.fasta
├── 4527_PUR2_codon_pos.tsv
└── ...  (222 pairs = 444 files total)
```

### Species shortlist and subsetted alignments

`shortlist.txt` contains 16 focal crop and model species:

| Species | Notes |
|---------|-------|
| *Arabidopsis thaliana* | |
| *Brassica napus* | |
| *Streptanthus carinatus* | proxy for *S. tortuosus* (not in dataset) |
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

For each of the 222 high-confidence genes, alignments were subsetted to these species only and copied alongside their annotation TSVs to `results/shortlist_files/`. 44 genes are missing at least one shortlist species in their alignment (most commonly *Streptanthus carinatus*, absent from 24 genes); those files contain fewer than 16 sequences. The annotation TSVs are copied as-is since column indices are independent of which species are present.

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
│   └── validate_translations.py
├── requirements.txt
└── README.md
```
