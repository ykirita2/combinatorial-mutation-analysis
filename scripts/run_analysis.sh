#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Combinatorial mutation library analysis
#
# Usage:
#   run_analysis.sh FASTQ REF_FASTA BED SPIKEIN_FASTA OUTPUT_PREFIX
#
# Example:
#   run_analysis.sh \
#       sample.fastq.gz \
#       target.fasta \
#       target.bed \
#       lambda_3.6kb.fasta \
#       sample1
# ============================================================

FASTQ=$1
REF=$2
BED=$3
SPIKEIN=$4
PREFIX=$5

THREADS=16
QC_THREADS=8
BGZIP_THREADS=4


# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

for f in "$FASTQ" "$REF" "$BED" "$SPIKEIN"
do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: File not found: $f" >&2
        exit 1
    fi
done


# ------------------------------------------------------------
# Reference information
# ------------------------------------------------------------

# Create FASTA index if absent
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

# This pipeline assumes a single-sequence amplicon FASTA
CONTIG=$(cut -f1 "${REF}.fai" | head -n 1)
REF_LENGTH=$(cut -f2 "${REF}.fai" | head -n 1)

echo "============================================"
echo "Combinatorial mutation library analysis"
echo "============================================"
echo "FASTQ      : $FASTQ"
echo "Reference  : $REF"
echo "Contig     : $CONTIG"
echo "Length     : $REF_LENGTH"
echo "BED        : $BED"
echo "Output     : $PREFIX"
echo "Start      : $(date)"
echo "============================================"


# ------------------------------------------------------------
# (1) Raw read statistics
# ------------------------------------------------------------

NanoPlot \
    --drop_outliers \
    -t "$QC_THREADS" \
    -o "${PREFIX}.NanoPlot" \
    --fastq "$FASTQ"


# ------------------------------------------------------------
# (2) Read QC
# Low-quality read and spike-in filtering
# ------------------------------------------------------------

zcat "$FASTQ" |
    chopper \
        -q 10 \
        -t "$QC_THREADS" \
        -c "$SPIKEIN" |
    bgzip \
        -c \
        -@ "$BGZIP_THREADS" \
    > "${PREFIX}.clean.fastq.gz"


# ------------------------------------------------------------
# (3) Clean read statistics
# ------------------------------------------------------------

NanoPlot \
    --drop_outliers \
    -t "$THREADS" \
    -o "${PREFIX}.clean.NanoPlot" \
    --fastq "${PREFIX}.clean.fastq.gz"


# ------------------------------------------------------------
# (4) Mapping
# ------------------------------------------------------------

minimap2 \
    -t "$THREADS" \
    -ax map-hifi \
    "$REF" \
    "${PREFIX}.clean.fastq.gz" |
    samtools sort \
        -@ "$THREADS" \
        -o "${PREFIX}.clean.bam" \
        --reference "$REF" \
        -

samtools index "${PREFIX}.clean.bam"


# ------------------------------------------------------------
# (5) Full-length primary alignment extraction
# ------------------------------------------------------------

# Reads overlapping the first base
samtools view \
    -b \
    "${PREFIX}.clean.bam" \
    "${CONTIG}:1-1" \
    > "${PREFIX}.clean.left.bam"

samtools index "${PREFIX}.clean.left.bam"


# Reads also overlapping the last base
# -F 2304 removes secondary (256) + supplementary (2048)
samtools view \
    -F 2304 \
    -b \
    "${PREFIX}.clean.left.bam" \
    "${CONTIG}:${REF_LENGTH}-${REF_LENGTH}" \
    > "${PREFIX}.clean.full.bam"

samtools index "${PREFIX}.clean.full.bam"


# ------------------------------------------------------------
# (6) Sequence-pattern counting
# ------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

perl "${SCRIPT_DIR}/read_count_by_seq_pattern.pl" \
    "$REF" \
    "$BED" \
    "${PREFIX}.clean.full.bam" \
    > "${PREFIX}.out"


echo "============================================"
echo "Finished"
echo "End: $(date)"
echo "Output: ${PREFIX}.out"
echo "============================================"
