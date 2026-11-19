# Combinatorial Mutation Library Analysis

Nanopore sequencing pipeline for combinatorial mutation library analysis.

This repository contains the analysis scripts, software environment, and lambda spike-in reference sequence.  
Sample-specific FASTQ files, amplicon reference sequences, BED files, and analysis outputs are not stored in this repository.

## Repository structure

```text
nanopore-combinatorial-mutation-analysis/
├── README.md
├── environment.yml
├── .gitignore
├── reference/
│   └── lambda_3.6kb.fasta
└── scripts/
    ├── run_analysis.sh
    └── read_count_by_seq_pattern.pl
```

## Analysis workflow

The pipeline performs the following steps:

1. Raw Nanopore read QC using NanoPlot
2. Filtering of low-quality reads using Chopper
3. Removal of lambda spike-in reads
4. QC of filtered reads using NanoPlot
5. Mapping to the amplicon reference using minimap2
6. Extraction of reads spanning the full reference sequence
7. Removal of secondary and supplementary alignments
8. Counting reads according to sequence patterns specified in a BED file

## Required input files

The following files are supplied separately for each analysis and are not stored in this repository.

### Nanopore FASTQ

```text
sample.fastq.gz
```

### Amplicon reference

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

### Sequence-pattern BED file

```text
pattern.bed
```

The BED file specifies the regions and sequence patterns used for read counting.

Coordinates are **0-based**.

### Lambda spike-in reference

```text
reference/lambda_3.6kb.fasta
```

This file is included in the repository and is used by Chopper to remove spike-in control reads.

## Software environment

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

## Usage

Run the pipeline as:

```bash
bash /path/to/nanopore-combinatorial-mutation-analysis/scripts/run_analysis.sh \
    sample.fastq.gz \
    ref.fasta \
    pattern.bed \
    output_prefix
```

For example:

```bash
bash ~/github/nanopore-combinatorial-mutation-analysis/scripts/run_analysis.sh \
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

## Output

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

## Custom Perl script

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

## Data policy

Raw sequencing data and analysis-specific files are intentionally excluded from this repository.

The following files are not tracked by Git:

- FASTQ files
- Analysis-specific FASTA files
- BED files
- BAM/SAM files
- NanoPlot outputs
- Other intermediate analysis files

Only the common lambda spike-in reference is included in the repository.

## Notes

This pipeline assumes that the reference FASTA contains a single amplicon sequence.

Reads overlapping both the first and last bases of the reference are retained as full-length aligned reads. Secondary and supplementary alignments are excluded before sequence-pattern counting.

## Combinatorial epistasis analysis of Cas12j

The complete combinatorial epistasis analysis of Cas12j,
including genotype-level measurements, Python analysis
scripts, and derived results, is available in
[epistasis/](epistasis/README.md).

