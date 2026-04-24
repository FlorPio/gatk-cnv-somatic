# GATK CNV Somatic Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.04.0-23aa62.svg)](https://www.nextflow.io/)
[![nf-core](https://img.shields.io/badge/nf--core-template-23aa62.svg)](https://nf-co.re/)
[![Docker](https://img.shields.io/badge/docker-florpio%2Fcnv--references-blue.svg)](https://hub.docker.com/r/florpio/cnv-references)

A Nextflow pipeline for somatic Copy Number Variation (CNV) analysis in tumor-normal paired samples using GATK best practices.

## Overview

This pipeline performs CNV calling on tumor-normal matched pairs using the GATK4 somatic CNV workflow. It includes:

- Automatic interval preprocessing (accepts `.bed` or `.interval_list`)
- Read count collection and denoising with Panel of Normals (tumor only)
- Allelic count collection for tumor and matched normal
- Segmentation with allelic fraction modeling (ModelSegments)
- CNV calling with copy ratio thresholds
- Diagnostic plots with configurable minimum contig length
- Gene annotation with MANE transcripts and external gene list

## Pipeline Steps

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        GATK CNV Somatic Pipeline                        │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  Tumor BAM ──┐    Normal BAM ──┐    References:                        │
│              │                 │    • Panel of Normals (.hdf5)          │
│              │                 │    • Common SNPs (.vcf.gz)             │
│              │                 │    • Intervals (.bed or .interval_list)│
│              │                 │    • FASTA (main chroms only)          │
│              ▼                 │                                        │
│  ┌──────────────────────────┐  │                                        │
│  │ 0. PreprocessIntervals   │  │  ◄── only if .bed provided            │
│  │    (auto-detected)       │  │      skipped for .interval_list        │
│  └────────────┬─────────────┘  │                                        │
│               ▼                │                                        │
│  ┌──────────────────────────┐  │                                        │
│  │ 1. CollectReadCounts     │  │  ◄── tumor only                       │
│  └────────────┬─────────────┘  │                                        │
│               ▼                │                                        │
│  ┌──────────────────────────┐  │                                        │
│  │ 2. DenoiseReadCounts     │  │  ◄── tumor only, uses PON             │
│  └────────────┬─────────────┘  │                                        │
│               │                │                                        │
│               ▼                ▼                                        │
│  ┌─────────────────┐ ┌─────────────────┐                               │
│  │ 3. Collect      │ │ 3. Collect      │                               │
│  │ AllelicCounts   │ │ AllelicCounts   │ ◄── Common SNPs + FASTA       │
│  │ (Tumor)         │ │ (Normal)        │                               │
│  └────────┬────────┘ └────────┬────────┘                               │
│           │                   │                                        │
│           └─────────┬─────────┘                                        │
│                     ▼                                                   │
│  ┌──────────────────────────────────────┐                              │
│  │ 4. ModelSegments                      │                              │
│  │    Inputs: DenoisedCR(T) +            │                              │
│  │    AllelicCounts(T) + AllelicCounts(N)│                              │
│  │    Outputs: .modelFinal.seg           │                              │
│  │             .cr.seg                   │                              │
│  │             .hets.tsv                 │                              │
│  └──────────────────┬───────────────────┘                              │
│                     │                                                   │
│         ┌───────────┼──────────────┐                                   │
│         ▼           ▼              ▼                                   │
│  ┌───────────┐ ┌──────────┐ ┌──────────────┐                          │
│  │ 5. Call   │ │ 6. Plot  │ │ 7. Annotate  │                          │
│  │ CNVs      │ │ Results  │ │ Genes        │                          │
│  │ (.cr.seg) │ │ (+DICT)  │ │ (MANE+Genes) │                          │
│  └───────────┘ └──────────┘ └──────────────┘                          │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

> **Note:** Steps 1-2 run on tumor samples only. Normal BAMs are used only for CollectAllelicCounts (step 3).
> Normal QC (denoising, segmentation) is handled by the companion pipeline [gatk-cnv-pon](https://github.com/FlorPio/gatk-cnv-pon).

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 25.04.0
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)

## Quick Start

### 1. Setup References

Prepare all reference files with the setup script (run once per genome/capture kit):

```bash
./docs/setup_references.sh \
    --fasta /path/to/hg38.fa \
    --bed /path/to/capture_targets.bed \
    --output-dir references
```

This creates the filtered FASTA (main chromosomes only), preprocessed intervals, annotated intervals (GC content for PON), and common SNPs from 1000 Genomes.

Alternatively, the reference files are available in the Docker image:

```bash
docker pull florpio/cnv-references:hg38-v1.0
docker create --name cnv-refs-temp florpio/cnv-references:hg38-v1.0
docker cp cnv-refs-temp:/references ./references
docker rm cnv-refs-temp
```

### 2. Prepare Samplesheet

Create a CSV file with your tumor-normal pairs:

```csv
sample,tumor_bam,normal_bam
Patient_01,/path/to/Patient_01_T01.recal.bam,/path/to/Patient_01_N01.recal.bam
Patient_02,/path/to/Patient_02_T01.recal.bam,/path/to/Patient_02_N01.recal.bam
```

**Requirements:**
- BAM files must be aligned to hg38 (main chromosomes only)
- BAM files must be sorted and indexed (.bai)
- Tumor and normal must be from the same patient

### 3. Create params.json

```json
{
  "fasta": "/path/to/references/genome/hg38.main_chroms.fa",
  "fai": "/path/to/references/genome/hg38.main_chroms.fa.fai",
  "dict": "/path/to/references/genome/hg38.main_chroms.dict",
  "intervals": "/path/to/references/intervals/preprocessed.interval_list",
  "pon": "/path/to/references/pon/wes-do-gc.pon.hdf5",
  "common_snps": "/path/to/references/snps/common_snps_preprocessed.vcf.gz",
  "mane_annotation": "/path/to/references/annotation/exons_mane.txt"
}
```

> **Note:** `--intervals` accepts both `.interval_list` (used directly) and `.bed` files (auto-triggers `PreprocessIntervals`). If you provide a `.bed` file, ensure the same padding was used for PON creation.

### 4. Run Pipeline

```bash
nextflow run FlorPio/gatk-cnv-somatic \
    -profile docker \
    -params-file params.json \
    --input samplesheet.csv \
    --outdir results
```

## Intervals: `.bed` vs `.interval_list`

The pipeline auto-detects the input type:

| Input | Behavior |
|-------|----------|
| `--intervals targets.interval_list` | Used directly. PreprocessIntervals is skipped. |
| `--intervals targets.bed` | PreprocessIntervals runs automatically with `--padding` and `--bin_length`. |

**Important:** The same intervals, padding, and bin_length must be used for both PON creation and case analysis. If your PON was built with `--padding 50`, pass `--padding 50` when using a `.bed` file.

## Gene Annotation

CNV segments are annotated with cancer-relevant genes using MANE Select transcripts. The gene list is externalized to `assets/cancer_genes.txt` (one gene per line), which can be customized:

```bash
# Use the default list (~107 cancer genes)
nextflow run FlorPio/gatk-cnv-somatic ...

# Use a custom gene list
nextflow run FlorPio/gatk-cnv-somatic --genes_list /path/to/my_genes.txt ...
```

The annotation R script (`bin/annotate_cnvs.R`) includes tumor purity estimation via LOH-based inference and absolute copy number calculation.

## Output

```
results/
├── intervals/
│   └── preprocessed.interval_list       # Only if .bed input was used
├── counts/
│   ├── {tumor}_tumor.hdf5               # Raw read counts (tumor only)
│   ├── {tumor}_standardizedCR.tsv       # Standardized copy ratios
│   ├── {tumor}_denoisedCR.tsv           # Denoised copy ratios
│   ├── {tumor}.allelicCounts.tsv        # Tumor allelic counts
│   └── {normal}.allelicCounts.tsv       # Normal allelic counts
├── model_segments/
│   ├── {tumor}.modelFinal.seg           # Final segmentation
│   ├── {tumor}.cr.seg                   # Copy ratio segments
│   └── {tumor}.hets.tsv                 # Heterozygous sites
├── called_segments/
│   └── {tumor}.called.seg              # Called CNVs (+, -, 0)
├── plots/
│   ├── {tumor}.denoised.png            # Copy ratio plots
│   └── {tumor}.modeled.png             # Segmentation plots
├── annotation/
│   ├── {tumor}_annotated.txt           # Gene-annotated CNVs
│   └── {tumor}_purity.txt             # Tumor purity estimates
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    └── software_versions.yml
```

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--input` | Samplesheet CSV (sample, tumor_bam, normal_bam) |
| `--outdir` | Output directory |
| `--fasta` | Reference genome FASTA (main chromosomes only) |
| `--intervals` | Target intervals (`.bed` or `.interval_list`) |
| `--pon` | Panel of Normals HDF5 file |
| `--common_snps` | Common SNPs VCF (.vcf.gz with .tbi index) |

### Annotation

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--mane_annotation` | MANE Select exon annotation file | Required if annotation enabled |
| `--genes_list` | Cancer gene list (one gene per line) | `assets/cancer_genes.txt` |

### Interval preprocessing

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--padding` | Padding in bp (only used with `.bed` input) | `50` |
| `--bin_length` | Bin length (0 = no binning for WES) | `0` |
| `--exclude_intervals` | Blacklist regions to exclude | `null` |

### Plot options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--somatic_min_contig_length` | Min contig length for plots (excludes chrM) | `46709983` (chr22 hg38) |

### Skip options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--skip_plots` | Skip diagnostic plots | `false` |
| `--skip_annotation` | Skip gene annotation | `false` |

## Companion Pipelines

| Pipeline | Purpose |
|----------|---------|
| [gatk-cnv-pon](https://github.com/FlorPio/gatk-cnv-pon) | Create and manage the Panel of Normals |
| `docs/setup_references.sh` | One-time reference file preparation |

## Reference Files

**IMPORTANT:** The FASTA reference must contain main chromosomes only (chr1-22, chrX, chrY) without alt contigs, decoys, HLA, or unplaced scaffolds. Use `docs/setup_references.sh` to create a filtered reference from a full hg38 FASTA.

| File | Description | Created by |
|------|-------------|------------|
| hg38.main_chroms.fa | Reference genome (main chromosomes) | `setup_references.sh` |
| preprocessed.interval_list | Padded target intervals | `setup_references.sh` |
| annotated_intervals.tsv | GC content per interval (for PON) | `setup_references.sh` |
| common_snps_preprocessed.vcf.gz | 1000G SNPs filtered to targets | `setup_references.sh` |
| wes-do-gc.pon.hdf5 | Panel of Normals | `gatk-cnv-pon` pipeline |
| exons_mane.txt | MANE transcript exon coordinates | Provided |

## Troubleshooting

### Common Issues

**Error: Reference file does not exist**
- Run `docs/setup_references.sh` to prepare reference files
- Verify paths in `params.json` are absolute paths

**Error: Interval list does not match PON**
- Ensure the same intervals, padding, and bin_length were used for PON creation and case analysis

**Error: Out of memory**
- Increase memory in `conf/base.config` or use `-c custom.config`
- Use `-resume` to continue from last checkpoint

**Error: BAM index not found**
- Ensure `.bai` files exist alongside BAM files
- Run `samtools index` if missing

**Annotation failed but pipeline continued**
- This is expected behavior (fallback). Check logs for R script errors
- An empty annotated file with headers is created so downstream steps are not affected

### Resume Failed Run

```bash
nextflow run . -profile docker -params-file params.json --input samplesheet.csv --outdir results -resume
```

## Citation

If you use this pipeline, please cite:

- [GATK](https://gatk.broadinstitute.org/) - McKenna et al., 2010
- [Nextflow](https://www.nextflow.io/) - Di Tommaso et al., 2017
- [nf-core](https://nf-co.re/) - Ewels et al., 2020

## Author

**Florencia Piovaroli**
Bioinformatics | Cancer Genomics | CNV Analysis

- GitHub: [@FlorPio](https://github.com/FlorPio)

## License

MIT License - see [LICENSE](LICENSE) for details.
