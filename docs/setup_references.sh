#!/usr/bin/env bash
# ============================================================
# setup_references.sh
# ============================================================
# One-time reference preparation for the GATK CNV Somatic pipeline.
# Creates all required reference files from a full hg38 FASTA
# and a capture kit BED file.
#
# This script is meant to be run ONCE when setting up the
# analysis environment. The output files are reused across
# all subsequent pipeline runs.
#
# Requirements:
#   - Docker (or local GATK 4.6.2.0 + samtools + wget)
#   - ~15 GB disk space for intermediate files
#   - Internet access (to download 1000G SNPs)
#
# Usage:
#   ./setup_references.sh \
#       --fasta /path/to/hg38.fa \
#       --bed /path/to/capture_targets.bed \
#       --output-dir /path/to/references
#
# Output structure:
#   references/
#   ├── genome/
#   │   ├── hg38.main_chroms.fa
#   │   ├── hg38.main_chroms.fa.fai
#   │   └── hg38.main_chroms.dict
#   ├── intervals/
#   │   ├── exomas.interval_list
#   │   ├── preprocessed.interval_list
#   │   └── annotated_intervals.tsv      ← used for GC correction in PON
#   └── snps/
#       ├── common_snps_preprocessed.vcf.gz
#       └── common_snps_preprocessed.vcf.gz.tbi
#
# NOTE: The Panel of Normals (PON) is NOT created by this
# script. Use the gatk-cnv-pon pipeline for that.
# When creating the PON, pass annotated_intervals.tsv via
# --annotated-intervals to enable explicit GC correction.
# ============================================================

set -euo pipefail

# ============================================================
# CONFIGURATION
# ============================================================

FASTA=""
BED=""
OUTPUT_DIR="references"
PADDING=50           # Must match PON creation padding
BIN_LENGTH=0         # 0 = no binning (WES). 1000 = WGS
GATK_IMAGE="broadinstitute/gatk:4.6.2.0"
SNPS_URL="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

# Main chromosomes to keep (no alt, decoy, HLA, unplaced)
MAIN_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

# ============================================================
# PARSE ARGUMENTS
# ============================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta)       FASTA="$2";       shift 2 ;;
        --bed)         BED="$2";         shift 2 ;;
        --output-dir)  OUTPUT_DIR="$2";  shift 2 ;;
        --padding)     PADDING="$2";     shift 2 ;;
        --help|-h)
            head -35 "$0" | tail -30
            return 0 2>/dev/null || true
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            return 1 2>/dev/null || true
            ;;
    esac
done

if [[ -z "$FASTA" ]] || [[ -z "$BED" ]]; then
    echo "ERROR: --fasta and --bed are required"
    echo "Run with --help for usage"
    return 1 2>/dev/null || true
fi

# ============================================================
# SETUP
# ============================================================

GENOME_DIR="${OUTPUT_DIR}/genome"
INTERVALS_DIR="${OUTPUT_DIR}/intervals"
SNPS_DIR="${OUTPUT_DIR}/snps"

mkdir -p "$GENOME_DIR" "$INTERVALS_DIR" "$SNPS_DIR"

# Docker wrapper
docker_gatk() {
    docker run --rm \
        -v "$(realpath "$OUTPUT_DIR"):$(realpath "$OUTPUT_DIR")" \
        -v "$(dirname "$(realpath "$FASTA")"):$(dirname "$(realpath "$FASTA")"):ro" \
        -v "$(dirname "$(realpath "$BED")"):$(dirname "$(realpath "$BED")"):ro" \
        -w "$(realpath "$OUTPUT_DIR")" \
        "${GATK_IMAGE}" gatk "$@"
}

echo ""
echo "╔══════════════════════════════════════════════╗"
echo "║   GATK CNV Somatic — Reference Setup        ║"
echo "╚══════════════════════════════════════════════╝"
echo ""
echo "  FASTA:      $FASTA"
echo "  BED:        $BED"
echo "  Output:     $OUTPUT_DIR"
echo "  Padding:    $PADDING bp"
echo ""

# ============================================================
# STEP 1: Create main-chromosomes-only reference
# ============================================================
# The GATK CNV workflow requires a reference with only main
# chromosomes. Alt contigs, decoys, HLA, and unplaced scaffolds
# cause issues with the Panel of Normals and interval matching.
# ============================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 1/5: Creating main-chromosomes-only reference"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

MAIN_FA="${GENOME_DIR}/hg38.main_chroms.fa"

if [[ -f "$MAIN_FA" ]]; then
    echo "  Skipping (already exists)"
else
    samtools faidx "$FASTA" $MAIN_CHROMS > "$MAIN_FA"
    echo "  Created: $MAIN_FA"
fi

# ============================================================
# STEP 2: Index and create sequence dictionary
# ============================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 2/5: Creating FASTA index and sequence dictionary"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [[ -f "${MAIN_FA}.fai" ]]; then
    echo "  Index: already exists"
else
    samtools faidx "$MAIN_FA"
    echo "  Created: ${MAIN_FA}.fai"
fi

MAIN_DICT="${GENOME_DIR}/hg38.main_chroms.dict"
if [[ -f "$MAIN_DICT" ]]; then
    echo "  Dictionary: already exists"
else
    docker_gatk CreateSequenceDictionary \
        -R "$MAIN_FA" \
        -O "$MAIN_DICT"
    echo "  Created: $MAIN_DICT"
fi

# ============================================================
# STEP 3: BedToIntervalList + PreprocessIntervals
# ============================================================
# Converts the capture kit BED to GATK interval list format,
# then pads and processes intervals for CNV analysis.
#
# IMPORTANT: The same padding and bin_length values MUST be
# used when creating the Panel of Normals. Record these values.
# ============================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 3/6: Converting BED → interval list → preprocessed"
echo "  Padding: ${PADDING} bp | Bin length: ${BIN_LENGTH}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

INTERVAL_LIST="${INTERVALS_DIR}/exomas.interval_list"
PREPROCESSED="${INTERVALS_DIR}/preprocessed.interval_list"

if [[ -f "$INTERVAL_LIST" ]]; then
    echo "  BedToIntervalList: already exists"
else
    docker_gatk BedToIntervalList \
        -I "$BED" \
        -O "$INTERVAL_LIST" \
        -SD "$MAIN_DICT"
    echo "  Created: $INTERVAL_LIST"
fi

if [[ -f "$PREPROCESSED" ]]; then
    echo "  PreprocessIntervals: already exists"
else
    docker_gatk PreprocessIntervals \
        -L "$INTERVAL_LIST" \
        -R "$MAIN_FA" \
        --padding "$PADDING" \
        --bin-length "$BIN_LENGTH" \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O "$PREPROCESSED"
    echo "  Created: $PREPROCESSED"
fi

# ============================================================
# STEP 4: AnnotateIntervals (GC content for PON correction)
# ============================================================
# Computes GC content per interval. This is used by
# CreateReadCountPanelOfNormals via --annotated-intervals
# to enable explicit GC correction in the PON.
#
# The GC info is baked into the PON, so DenoiseReadCounts
# does NOT need to be passed annotated_intervals.tsv when
# using a PON (would cause double correction).
#
# Optional but recommended: improves denoising quality,
# especially for capture kits with variable GC bias.
# ============================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 4/6: Annotating intervals with GC content"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

ANNOTATED="${INTERVALS_DIR}/annotated_intervals.tsv"

if [[ -f "$ANNOTATED" ]]; then
    echo "  Skipping (already exists)"
else
    docker_gatk AnnotateIntervals \
        -L "$PREPROCESSED" \
        -R "$MAIN_FA" \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O "$ANNOTATED"
    echo "  Created: $ANNOTATED"
fi

# ============================================================
# STEP 5: Download common SNPs (1000 Genomes)
# ============================================================
# Source: 1000 Genomes Phase 1 high-confidence SNPs (hg38)
# These are used by CollectAllelicCounts to determine
# allelic fractions at known heterozygous sites.
# ============================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 5/6: Downloading common SNPs from 1000 Genomes"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

RAW_SNPS="${SNPS_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

if [[ -f "$RAW_SNPS" ]]; then
    echo "  Skipping download (already exists)"
else
    wget -q --show-progress -O "$RAW_SNPS" "$SNPS_URL"
    wget -q --show-progress -O "${RAW_SNPS}.tbi" "${SNPS_URL}.tbi"
    echo "  Downloaded: $RAW_SNPS"
fi

# ============================================================
# STEP 5: Filter SNPs to target intervals
# ============================================================
# Restricts the common SNPs to only those within the
# preprocessed target intervals. This reduces the file
# from ~7M sites to ~100-300K (for typical WES kits).
# ============================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 6/6: Filtering SNPs to target intervals"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

FILTERED_SNPS="${SNPS_DIR}/common_snps_preprocessed.vcf.gz"

if [[ -f "$FILTERED_SNPS" ]]; then
    echo "  Skipping (already exists)"
else
    docker_gatk SelectVariants \
        -V "$RAW_SNPS" \
        -L "$PREPROCESSED" \
        -O "$FILTERED_SNPS"
    echo "  Created: $FILTERED_SNPS"
fi

# ============================================================
# SUMMARY
# ============================================================

echo ""
echo "╔══════════════════════════════════════════════╗"
echo "║   Reference setup complete                   ║"
echo "╚══════════════════════════════════════════════╝"
echo ""
echo "  Output files:"
echo "    ${MAIN_FA}"
echo "    ${MAIN_FA}.fai"
echo "    ${MAIN_DICT}"
echo "    ${PREPROCESSED}"
echo "    ${ANNOTATED}"
echo "    ${FILTERED_SNPS}"
echo "    ${FILTERED_SNPS}.tbi"
echo ""
echo "  Use in params.json (gatk-cnv-somatic):"
echo "    {"
echo "      \"fasta\":       \"${MAIN_FA}\","
echo "      \"fai\":         \"${MAIN_FA}.fai\","
echo "      \"dict\":        \"${MAIN_DICT}\","
echo "      \"intervals\":   \"${PREPROCESSED}\","
echo "      \"common_snps\": \"${FILTERED_SNPS}\","
echo "      \"pon\":         \"/path/to/your/pon.hdf5\""
echo "    }"
echo ""
echo "  IMPORTANT: Use the SAME padding (${PADDING}) when"
echo "  creating the Panel of Normals with gatk-cnv-pon."
echo ""
echo "  When creating the PON, pass annotated_intervals.tsv"
echo "  to enable explicit GC correction:"
echo "    --annotated_intervals ${ANNOTATED}"
echo ""
echo "  Note: do NOT pass annotated_intervals to the somatic"
echo "  pipeline — GC correction is already baked into the PON."
echo ""
echo "  Next step: Create a PON with gatk-cnv-pon pipeline"
echo "  or use pon_manager.sh to manage normal samples."
echo ""
