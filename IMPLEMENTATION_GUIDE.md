# gatk-cnv-somatic — Implementation Guide v1.2.0

## Quick Overview

You have 4 types of changes:
- **NEW files** → create them
- **REPLACE files** → overwrite entirely
- **DELETE references** → remove parameter that no longer exists
- **INSTALL modules** → nf-core modules already in your repo (no change needed)

---

## Step 1: Clone and create branch

```bash
git clone https://github.com/FlorPio/gatk-cnv-somatic.git
cd gatk-cnv-somatic
git checkout -b feature/v1.2.0-optimization
```

---

## Step 2: NEW files to create

### 2a. `assets/cancer_genes.txt`
Cancer gene list externalized from R script.
```bash
# File provided in outputs: assets/cancer_genes.txt
cp /path/to/downloaded/assets/cancer_genes.txt assets/cancer_genes.txt
git add assets/cancer_genes.txt
```

### 2b. `bin/annotate_cnvs.R`
R script moved from Docker image to pipeline repo.
```bash
mkdir -p bin
cp /path/to/downloaded/bin/annotate_cnvs.R bin/annotate_cnvs.R
chmod +x bin/annotate_cnvs.R
git add bin/annotate_cnvs.R
```

### 2c. `CHANGELOG.md`
```bash
cp /path/to/downloaded/CHANGELOG.md CHANGELOG.md
git add CHANGELOG.md
```

---

## Step 3: REPLACE files (overwrite entirely)

### 3a. `workflows/gatk-cnv-somatic.nf` ← MAIN CHANGES HERE
Changes:
- CollectReadCounts/DenoiseReadCounts now run TUMOR ONLY
- Normal BAM only goes to CollectAllelicCounts
- Double-key joins (meta.id for tumor-tumor, combine for tumor-normal)
- Annotation uses genes_list instead of annotate_script
- Fallback annotation (pipeline doesn't crash if R fails)
- minimum-contig-length passed via ext.args

```bash
cp /path/to/downloaded/final-somatic/workflows/gatk-cnv-somatic.nf workflows/gatk-cnv-somatic.nf
```

### 3b. `modules/local/annotate_cnv.nf`
Changes:
- Removed path(rscript) input
- Added path(genes_list) input
- Added conda directive
- Fallback: creates empty output if R fails
- Fixed duplicate output channels

```bash
cp /path/to/downloaded/modules/local/annotate_cnv.nf modules/local/annotate_cnv.nf
```

### 3c. `conf/modules.config`
Changes:
- publishDir for all modules
- GATK 4.6.2.0 container overrides for local modules
- --minimum-contig-length for plot modules

```bash
cp /path/to/downloaded/final-somatic/conf/modules.config conf/modules.config
```

### 3d. `nextflow.config`
Changes:
- Removed params.annotate_script
- Added params.genes_list
- Added params.somatic_min_contig_length
- Added params.publish_dir_mode
- FASTA documented as "main chroms only"
- Execution reports enabled (timeline, report, trace, dag)

```bash
cp /path/to/downloaded/nextflow.config nextflow.config
```

### 3e. `nextflow_schema.json`
Changes:
- Removed annotate_script definition
- Added genes_list with help_text
- Added somatic_min_contig_length
- Improved fasta description

```bash
cp /path/to/downloaded/nextflow_schema.json nextflow_schema.json
```

### 3f. `params.template.json`
Changes:
- Removed annotate_script
- Added genes_list

```bash
cp /path/to/downloaded/params.template.json params.template.json
```

---

## Step 4: Files that DON'T change

These files stay exactly as they are:

| File | Why |
|---|---|
| `main.nf` (root) | Only orchestrates — doesn't call ANNOTATE_CNV |
| `subworkflows/local/main.nf` | Pipeline initialization, untouched |
| `subworkflows/nf-core/*` | nf-core boilerplate, untouched |
| `modules/nf-core/gatk4/collectreadcounts/` | nf-core module, untouched |
| `modules/nf-core/gatk4/denoisereadcounts/` | nf-core module, untouched |
| `modules/local/modelsegments_somatic.nf` | Container overridden via modules.config |
| `modules/local/callcopyratiosegments.nf` | Container overridden via modules.config |
| `modules/local/collectalleliccounts.nf` | Container overridden via modules.config |
| `modules/local/plotdenoisedcopyratios.nf` | Container overridden via modules.config |
| `modules/local/plotmodeledsegments.nf` | Container overridden via modules.config |
| `assets/schema_input.json` | Samplesheet schema unchanged |
| `.gitignore` | Unchanged |
| `.nf-core.yml` | Unchanged |
| `modules.json` | Unchanged |

Note: GATK version for local modules is overridden to 4.6.2.0
via conf/modules.config container directives — the .nf files
themselves keep their original container declaration as fallback.

---

## Step 5: Git commands

```bash
# Stage all changes
git add -A

# Review what changed
git status
git diff --cached --stat

# Commit
git commit -m "feat: optimize pipeline for nf-core compliance (v1.2.0)

- Move R annotation script to bin/ (nf-core standard)
- Externalize gene list to assets/cancer_genes.txt
- Add --genes_list parameter, remove --annotate_script
- Run CollectReadCounts/DenoiseReadCounts on tumor only
  (normal only needs CollectAllelicCounts for ModelSegments)
- Robust double-key joins for multi-tumor support
- Add annotation fallback (pipeline continues if R fails)
- Add --somatic_min_contig_length for diagnostic plots
- Centralize publishDir and GATK 4.6.2.0 in conf/modules.config
- Document FASTA as main-chromosomes-only requirement
- Add CHANGELOG.md"

# Push
git push origin feature/v1.2.0-optimization

# Create PR on GitHub or merge directly
# Option A: PR (recommended)
#   Go to GitHub → Pull Requests → New → feature/v1.2.0-optimization → main

# Option B: Direct merge
git checkout main
git merge feature/v1.2.0-optimization
git push origin main
git tag -a v1.2.0 -m "v1.2.0: nf-core optimization + external gene list"
git push origin v1.2.0
```

---

## Step 6: Post-merge tasks

### Update your local params.json
```json
{
    "input": "samplesheet.csv",
    "outdir": "results",
    "fasta": "/home/golu/supporting_files/cnv/main_chroms/hg38.main_chroms.fa",
    "fai": "/home/golu/supporting_files/cnv/main_chroms/hg38.main_chroms.fa.fai",
    "dict": "/home/golu/supporting_files/cnv/main_chroms/hg38.main_chroms.dict",
    "intervals": "/home/golu/supporting_files/cnv/preprocessed.interval_list",
    "pon": "/home/golu/supporting_files/cnv/wes-do-gc.pon.hdf5",
    "common_snps": "/home/golu/supporting_files/cnv/common_snps_preprocessed.vcf.gz",
    "mane_annotation": "/home/golu/supporting_files/cnv/exons_mane.txt"
}
```
Note: `annotate_script` REMOVED, `genes_list` uses default from assets/.

### Test run
```bash
nextflow run FlorPio/gatk-cnv-somatic \
    -profile docker \
    -params-file params.json \
    --input samplesheet.csv \
    --outdir results_v1.2.0_test \
    -resume
```

### Verify output structure
```
results_v1.2.0_test/
├── counts/
│   ├── {tumor}_standardizedCR.tsv
│   ├── {tumor}_denoisedCR.tsv
│   └── {tumor|normal}.allelicCounts.tsv
├── model_segments/
│   ├── {tumor}.modelFinal.seg
│   ├── {tumor}.cr.seg
│   └── {tumor}.hets.tsv
├── called_segments/
│   └── {tumor}.called.seg
├── plots/
│   ├── {tumor}.denoised.png
│   └── {tumor}.modeled.png
├── annotation/
│   ├── {tumor}_annotated.txt
│   └── {tumor}_purity.txt
└── pipeline_info/
    └── gatk_cnv_somatic_software_versions.yml
```

---

## Files provided in this session (download checklist)

| Downloaded file | Goes to | Action |
|---|---|---|
| `assets/cancer_genes.txt` | `assets/cancer_genes.txt` | NEW |
| `bin/annotate_cnvs.R` | `bin/annotate_cnvs.R` | NEW |
| `CHANGELOG.md` | `CHANGELOG.md` | NEW |
| `final-somatic/workflows/gatk-cnv-somatic.nf` | `workflows/gatk-cnv-somatic.nf` | REPLACE |
| `final-somatic/conf/modules.config` | `conf/modules.config` | REPLACE |
| `modules/local/annotate_cnv.nf` | `modules/local/annotate_cnv.nf` | REPLACE |
| `nextflow.config` | `nextflow.config` | REPLACE |
| `nextflow_schema.json` | `nextflow_schema.json` | REPLACE |
| `params.template.json` | `params.template.json` | REPLACE |
