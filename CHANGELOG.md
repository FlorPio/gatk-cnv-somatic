# FlorPio/gatk-cnv-somatic: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-02-23

### Breaking Changes

- **Removed `params.annotate_script`** — The R annotation script is no longer passed as a pipeline parameter. It now lives in `bin/annotate_cnvs.R` following nf-core conventions, and is automatically available in the process `$PATH`. Users who referenced `annotate_script` in their `params.json` must remove it.

### Added

- **External gene list file** (`assets/cancer_genes.txt`) — The list of ~107 cancer-relevant genes (tumor suppressors and oncogenes) previously hardcoded in the R script is now an external text file. Users can supply their own custom gene list via `--genes_list` to filter CNV annotation results without modifying pipeline code.
- **`--genes_list` parameter** — New pipeline parameter pointing to the gene list file. Defaults to `${projectDir}/assets/cancer_genes.txt`.
- **`conf/modules.config`** — Centralized `publishDir` configuration for all modules and GATK container version overrides, following nf-core standards.
- **`conda` directive in `ANNOTATE_CNV`** — Module now declares its conda dependencies (`r-base`, `r-optparse`, `r-dplyr`, `r-tidyr`, `bioconductor-genomicranges`), enabling execution without Docker.
- **Execution reports** — `timeline`, `report`, `trace`, and `dag` are now enabled by default in `nextflow.config`.
- **`params.publish_dir_mode`** — Configurable publish directory mode (default: `copy`).

### Changed

- **R annotation script moved to `bin/`** — `annotate_cnvs_complete_fixed.R` renamed to `bin/annotate_cnvs.R`. Nextflow automatically adds `bin/` to `$PATH` in all processes, eliminating the need to pass the script as an input file.
- **R script accepts `--genes_list` parameter** — Replaces the hardcoded `genes_list <- c(...)` block (previously lines 78-94). The script now reads gene symbols from an external file, one per line, with `#` comment support.
- **R script comments translated to English** — All Spanish comments and log messages converted to English for international visibility and consistency.
- **GATK version standardized to 4.6.2.0** — All local modules (`MODELSEGMENTS_SOMATIC`, `CALLCOPYRATIOSEGMENTS`, `COLLECTALLELICCOUNTS`, `PLOTDENOISEDCOPYRATIOS`, `PLOTMODELEDSEGMENTS`) updated from `gatk:4.4.0.0` to `gatk:4.6.2.0` via `conf/modules.config` container overrides. Version 4.6.2.0 was validated to produce improved segmentation accuracy.
- **`ANNOTATE_CNV` module refactored**:
  - Input `path(rscript)` replaced with `path(genes_list)`.
  - Fixed duplicate output channels: `annotated` and `summary` previously both emitted `*.txt`. Now `annotated` emits `*_annotated.txt` and `purity` emits `*_purity.txt`.
  - Output file renaming logic improved to handle R script outputs correctly.
- **`nextflow_schema.json` updated** — Added `annotation_options` group with `genes_list` parameter definition. Removed `annotate_script`. Added `help_text` descriptions for reference file parameters.
- **`params.template.json` updated** — Replaced `annotate_script` with `genes_list`.

### Fixed

- **Duplicate output channels in `ANNOTATE_CNV`** — The `annotated` and `summary` emit channels both matched `*.txt`, causing Nextflow to emit the same files on both channels. Now each channel has a distinct glob pattern.

### Removed

- **`params.annotate_script`** — No longer needed; the R script is bundled in `bin/`.
- **Hardcoded gene list in R script** — Replaced by external `assets/cancer_genes.txt`.

---

## [1.0.0] - 2025-XX-XX

### Added

- Initial release of the GATK CNV Somatic pipeline.
- Nextflow DSL2 implementation based on nf-core template.
- Support for tumor-normal paired CNV analysis using GATK4 best practices.
- Steps: CollectReadCounts, DenoiseReadCounts, CollectAllelicCounts, ModelSegments, CallCopyRatioSegments.
- Optional diagnostic plots (PlotDenoisedCopyRatios, PlotModeledSegments).
- Optional gene annotation with MANE Select transcripts and purity-corrected absolute copy number.
- Docker images for references (`florpio/cnv-references:hg38-v1.0`) and R annotation (`florpio/cnv-annotate-r:1.0`).
- Template configuration files and example samplesheet.
