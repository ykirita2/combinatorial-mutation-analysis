# Combinatorial Mutation Library Analysis

This repository contains code and supporting resources for two related components of combinatorial mutation-library analysis:

1. **Nanopore sequencing analysis:** quality control, amplicon alignment, selection of full-length reads, and counting of designed sequence patterns.
2. **Cas12j combinatorial epistasis analysis:** analysis of precomputed genotype-level measurements, quantification of mutation interactions, fitness-landscape structure, model-order comparison, and in silico engineering strategies.

The Nanopore workflow uses sample-specific inputs supplied by the user; raw FASTQ files and sample-specific intermediates are **not** included here. The Cas12j analysis includes curated genotype-level CSV measurements and derived result tables in [`epistasis/`](epistasis/README.md). The two workflows are provided as separate components; the epistasis script does not automatically run the Nanopore pipeline.

## Repository structure

```text
combinatorial-mutation-analysis/
├── README.md
├── environment.yml
├── .gitignore
├── reference/
│   └── lambda_3.6kb.fasta
├── scripts/
│   ├── run_analysis.sh
│   └── read_count_by_seq_pattern.pl
└── epistasis/
    ├── README.md
    ├── data/
    │   ├── combinatorial_library_N.csv
    │   ├── combinatorial_library_C.csv
    │   └── substitution_assignment.csv
    ├── src/
    │   ├── epistasis_analysis.py
    │   ├── build_fig3.py
    │   └── build_figS3.py
    └── results/
        ├── library_summary.csv
        ├── pairwise_interaction_detectability.csv
        ├── interaction_vs_background_load.csv
        ├── local_optima.csv
        ├── substitution_background_effects.csv
        ├── model_order_cv.csv
        ├── coefficient_transfer.csv
        ├── strategy_benchmark.csv
        └── strategy_stability.csv
```

## 1. Nanopore sequencing and sequence-pattern counting

This component provides a Nanopore sequencing pipeline for combinatorial mutation-library analysis. It contains the analysis scripts, software environment, and lambda spike-in reference sequence; sample FASTQ files, analysis-specific amplicon references, BED files, and output files are supplied or generated locally rather than tracked in Git.

### Analysis workflow

The pipeline performs the following steps:

1. Raw Nanopore read QC using NanoPlot
2. Filtering of low-quality reads using Chopper
3. Removal of lambda spike-in reads
4. QC of filtered reads using NanoPlot
5. Mapping to the amplicon reference using minimap2
6. Extraction of reads spanning the full reference sequence
7. Removal of secondary and supplementary alignments
8. Counting reads according to sequence patterns specified in a BED file

### Required input files

The following files are supplied separately for each analysis and are not stored in this repository.

#### Nanopore FASTQ

```text
sample.fastq.gz
```

#### Amplicon reference

```text
ref.fasta
```

The reference FASTA should contain a single amplicon sequence.

Example:

```text
>target
ACGT...
```

The reference name and sequence length are automatically obtained from the FASTA index.

#### Sequence-pattern BED file

```text
pattern.bed
```

The BED file specifies the regions and sequence patterns used for read counting.

Coordinates are **0-based**.

#### Lambda spike-in reference

```text
reference/lambda_3.6kb.fasta
```

This file is included in the repository and is used by Chopper to remove spike-in control reads.

### Software environment

The required software environment is stored in:

```text
environment.yml
```

Create the conda environment with:

```bash
conda env create -f environment.yml
```

Then activate the environment:

```bash
conda activate <environment_name>
```

Main dependencies include:

- NanoPlot
- Chopper
- minimap2
- samtools
- bgzip / htslib
- Perl

### Usage

Run the pipeline as:

```bash
bash /path/to/combinatorial-mutation-analysis/scripts/run_analysis.sh \
    sample.fastq.gz \
    ref.fasta \
    pattern.bed \
    output_prefix
```

For example:

```bash
bash ~/github/combinatorial-mutation-analysis/scripts/run_analysis.sh \
    669_H1.fastq.gz \
    cas9d.fasta \
    cas9d_v4.bed \
    669_H1
```

The lambda spike-in reference does not need to be specified because the pipeline automatically uses:

```text
reference/lambda_3.6kb.fasta
```

from this repository.

### Output

The pipeline generates files with the specified output prefix.

Example:

```text
669_H1.NanoPlot/
669_H1.clean.fastq.gz
669_H1.clean.NanoPlot/
669_H1.clean.bam
669_H1.clean.bam.bai
669_H1.clean.left.bam
669_H1.clean.left.bam.bai
669_H1.clean.full.bam
669_H1.clean.full.bam.bai
669_H1.out
```

The final sequence-pattern count is stored in:

```text
669_H1.out
```

### Custom Perl script

Sequence-pattern counting is performed using:

```text
scripts/read_count_by_seq_pattern.pl
```

The script takes three arguments:

```bash
read_count_by_seq_pattern.pl \
    ref.fasta \
    pattern.bed \
    full_length_reads.bam
```

### Data policy

For the Nanopore preprocessing workflow, raw sequencing data and sample-specific analysis files are intentionally excluded from this repository. This policy does **not** exclude the curated, genotype-level Cas12j measurements in `epistasis/data/` or the derived summary tables in `epistasis/results/`.

The following files are not tracked by Git:

- FASTQ files
- Analysis-specific FASTA files
- BED files
- BAM/SAM files
- NanoPlot outputs
- Other intermediate analysis files

For the Nanopore processing pipeline, only the common lambda spike-in reference is provided in the repository; sample-specific reference FASTA and BED files are not included.

### Notes

This pipeline assumes that the reference FASTA contains a single amplicon sequence.

Reads overlapping both the first and last bases of the reference are retained as full-length aligned reads. Secondary and supplementary alignments are excluded before sequence-pattern counting.

## 2. Combinatorial epistasis analysis of Cas12j

This component provides the genotype-level input tables, analysis script, and derived outputs for the combinatorial-library analyses described in *Complete combinatorial mutagenesis resolves antagonistic epistasis and yields a hypercompact Cas12j genome editor* (Toyoshima et al.). See the [dedicated epistasis README](epistasis/README.md) for definitions and analysis details.

### Experimental libraries and input data

Two combinatorial libraries are included:

| Library | Substitutions | Designed genotypes | Genotypes retained after read-count filtering |
| --- | ---: | ---: | ---: |
| N-terminal segment | 11 | 2,048 | 1,823 |
| C-terminal segment | 10 | 1,024 | 1,021 |

The input files in `epistasis/data/` contain the measured genotype scores and replicate measurements for the two libraries, and a mapping between design columns and amino-acid substitutions:

- `combinatorial_library_N.csv`
- `combinatorial_library_C.csv`
- `substitution_assignment.csv`

These are processed genotype-level data, **not** raw Nanopore FASTQ files. The analysis applies read-count thresholds of 888 (N segment) and 1,647 (C segment) as specified in the source script.

### Software requirements

The main epistasis analysis was written for **Python 3.13** and requires `numpy`, `pandas`, and `scipy`. These requirements are separate from the Conda environment documented above for the Nanopore pipeline.

For example, in a Python 3.13 environment:

```bash
python -m pip install numpy pandas scipy
```

The supplied figure scripts also import `matplotlib`, but require additional in-memory objects and styling helpers; they are not self-contained command-line programs.

### Run the epistasis analysis

From the repository root:

```bash
python epistasis/src/epistasis_analysis.py --data epistasis/data --out epistasis/results
```

Or from within `epistasis/`:

```bash
cd epistasis
python src/epistasis_analysis.py --data data --out results
```

The main analysis reads the three input CSV files and writes eight result tables. It evaluates model-free pairwise interaction contrasts, background dependence of substitution effects, local optima and attraction basins, model order by cross-validation, and simulated engineering strategies. Its random-number generators use fixed seeds for the cross-validation and strategy simulation.

### Output files

| File in `epistasis/results/` | Description |
| --- | --- |
| `library_summary.csv` | Read-count filtering, retained genotypes, replicate diagnostics, and score uncertainty estimates |
| `pairwise_interaction_detectability.csv` | Pairwise interaction contrasts on parental and other backgrounds, uncertainty and multiple-testing-adjusted statistics |
| `interaction_vs_background_load.csv` | Mean interaction contrast versus the number of additional substitutions in the background |
| `local_optima.csv` | Single-substitution local optima and their basins of attraction |
| `substitution_background_effects.csv` | Individual substitution effects across accessible backgrounds |
| `model_order_cv.csv` | Cross-validated performance of additive, pairwise, and third-order models |
| `coefficient_transfer.csv` | Comparison of pairwise coefficients from low-order and complete-library fits |
| `strategy_benchmark.csv` | In silico comparison of seven mutation-combination / selection strategies |
| `strategy_stability.csv` | Additional supplied strategy-simulation stability results for multiple seeds and run counts |

**Reproducibility note:** The first eight tables are generated by `epistasis/src/epistasis_analysis.py`. The accompanying `strategy_stability.csv` is a separately supplied result table and is **not** recreated by that script.

### Interpretation and figure-generation notes

A pairwise interaction is calculated using a four-genotype (double-mutant-cycle) contrast:

```text
epsilon = y_ij - y_i - y_j + y_0
```

Contrasts are evaluated only when all four relevant genotypes are available. The benchmark simulates measurement noise and a common remeasurement-based candidate-selection step; details are described in [`epistasis/README.md`](epistasis/README.md).

`epistasis/src/build_fig3.py` and `epistasis/src/build_figS3.py` contain code for the corresponding figure panels. **These files require notebook/session variables and helper functions described in their headers; running them directly as standalone scripts is not supported by the supplied files.** The numeric result tables can be regenerated independently with `epistasis_analysis.py`.

### Data availability

The accompanying epistasis documentation reports the underlying sequencing data under SRA BioProject **PRJNA1514873**, and structural deposits **EMD-83117** (cryo-EM map) and **PDB 44ZO** (model).

## Scope and data management

The root `scripts/` workflow and `epistasis/` analysis serve different stages of analysis. The published `epistasis/data/` files are curated genotype-level measurement tables; sample-level raw sequencing files and other large Nanopore intermediates remain excluded from this repository, as described above. Running the Cas12j analysis does not require the sample-specific input FASTQ, FASTA, or BED files used by the general Nanopore pipeline.
