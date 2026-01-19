# GATK CNV Somatic Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.04.0-23aa62.svg)](https://www.nextflow.io/)
[![nf-core](https://img.shields.io/badge/nf--core-template-23aa62.svg)](https://nf-co.re/)
[![Docker](https://img.shields.io/badge/docker-florpio%2Fcnv--references-blue.svg)](https://hub.docker.com/r/florpio/cnv-references)

A Nextflow pipeline for somatic Copy Number Variation (CNV) analysis in tumor-normal paired samples using GATK best practices.

## Overview

This pipeline performs CNV calling on tumor-normal matched pairs using the GATK4 somatic CNV workflow. It includes:

- Read count collection and denoising with Panel of Normals
- Allelic count collection for tumor and matched normal
- Segmentation with allelic fraction modeling (ModelSegments)
- CNV calling with copy ratio thresholds
- Diagnostic plots generation
- Gene annotation with MANE transcripts

## Pipeline Steps

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        GATK CNV Somatic Pipeline                        │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────────────────┐  │
│  │ Tumor BAM    │    │ Normal BAM   │    │ References               │  │
│  └──────┬───────┘    └──────┬───────┘    │ • Panel of Normals       │  │
│         │                   │            │ • Common SNPs            │  │
│         ▼                   ▼            │ • Intervals              │  │
│  ┌──────────────────────────────────┐    └──────────────────────────┘  │
│  │ 1. CollectReadCounts (T + N)     │                                  │
│  └──────────────┬───────────────────┘                                  │
│                 ▼                                                       │
│  ┌──────────────────────────────────┐                                  │
│  │ 2. DenoiseReadCounts (T + N)     │◄── Panel of Normals             │
│  └──────────────┬───────────────────┘                                  │
│                 │                                                       │
│         ┌───────┴───────┐                                              │
│         ▼               ▼                                              │
│  ┌─────────────┐ ┌─────────────┐                                       │
│  │ 3. Collect  │ │ 3. Collect  │                                       │
│  │ AllelicCnts │ │ AllelicCnts │◄── Common SNPs                       │
│  │   (Tumor)   │ │  (Normal)   │                                       │
│  └──────┬──────┘ └──────┬──────┘                                       │
│         │               │                                              │
│         └───────┬───────┘                                              │
│                 ▼                                                       │
│  ┌──────────────────────────────────┐                                  │
│  │ 4. ModelSegments (Tumor-Normal)  │                                  │
│  └──────────────┬───────────────────┘                                  │
│                 │                                                       │
│         ┌───────┼───────┐                                              │
│         ▼       ▼       ▼                                              │
│  ┌───────┐ ┌────────┐ ┌────────────┐                                   │
│  │ 5.Call│ │ 6.Plot │ │ 7.Annotate │                                   │
│  │  CNVs │ │ Results│ │   Genes    │◄── MANE annotation               │
│  └───────┘ └────────┘ └────────────┘                                   │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 25.04.0
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)

## Quick Start

### 1. Setup References

The reference files are packaged in a Docker image for reproducibility.

```bash
# Download and run the setup script
./setup_references.sh

# Or manually:
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
- BAM files must be aligned to hg38
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
  "mane_annotation": "/path/to/references/annotation/exons_mane.txt",
  "annotate_script": "/path/to/references/annotation/annotate_cnvs_complete_fixed.R"
}
```

### 4. Run Pipeline

```bash
nextflow run FlorPio/gatk-cnv-somatic \
    -profile docker \
    -params-file params.json \
    --input samplesheet.csv \
    --outdir results
```

## Reference Files

The Docker image `florpio/cnv-references:hg38-v1.0` contains:

| File | Description | Path |
|------|-------------|------|
| hg38.main_chroms.fa | Reference genome (main chromosomes) | `/references/genome/` |
| preprocessed.interval_list | Target intervals for exome | `/references/intervals/` |
| wes-do-gc.pon.hdf5 | Panel of Normals for denoising | `/references/pon/` |
| common_snps_preprocessed.vcf.gz | Common SNPs for allelic counts | `/references/snps/` |
| exons_mane.txt | MANE transcript exon coordinates | `/references/annotation/` |
| annotate_cnvs_complete_fixed.R | CNV annotation script | `/references/annotation/` |

## Output

```
results/
├── counts/
│   ├── {sample}.hdf5                    # Raw read counts
│   ├── {sample}_standardizedCR.tsv      # Standardized copy ratios
│   └── {sample}_denoisedCR.tsv          # Denoised copy ratios
├── model_segments/
│   ├── {sample}.modelFinal.seg          # Final segmentation
│   ├── {sample}.cr.seg                  # Copy ratio segments
│   └── {sample}.hets.tsv                # Heterozygous sites
├── called_segments/
│   └── {sample}.called.seg              # Called CNVs (+, -, 0)
├── plots/
│   ├── {sample}.denoised.png            # Copy ratio plots
│   └── {sample}.modeled.png             # Segmentation plots
├── annotation/
│   └── {sample}_annotated.tsv           # Gene-annotated CNVs
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    └── software_versions.yml
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Samplesheet CSV | Required |
| `--outdir` | Output directory | Required |
| `--fasta` | Reference genome FASTA | Required |
| `--intervals` | Target intervals | Required |
| `--pon` | Panel of Normals | Required |
| `--common_snps` | Common SNPs VCF | Required |
| `--skip_plots` | Skip plot generation | `false` |
| `--skip_annotation` | Skip gene annotation | `false` |

## Pipeline Options

### Skip Steps

```bash
# Skip diagnostic plots
nextflow run . -profile docker --skip_plots true ...

# Skip gene annotation
nextflow run . -profile docker --skip_annotation true ...
```

### Resource Configuration

Edit `conf/base.config` or use `-c custom.config`:

```groovy
process {
    withName: 'MODELSEGMENTS_SOMATIC' {
        memory = '8.GB'
        cpus = 2
    }
}
```

## Troubleshooting

### Common Issues

**Error: Reference file does not exist**
- Run `./setup_references.sh` to extract reference files
- Verify paths in `params.json` are absolute paths

**Error: Out of memory**
- Reduce `params.max_memory` in `nextflow.config`
- Use `-resume` to continue from last checkpoint

**Error: BAM index not found**
- Ensure `.bai` files exist alongside BAM files
- Run `samtools index` if missing

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
