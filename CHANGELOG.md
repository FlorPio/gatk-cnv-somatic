# FlorPio/gatk-cnv-somatic: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] - 2026-04-09

### Breaking Changes

- **Removed `params.annotate_script`** — The R annotation script now lives in `bin/annotate_cnvs.R` following nf-core conventions and is automatically available in the process `$PATH`. Users who referenced `annotate_script` in their `params.json` must remove it.
- **CollectReadCounts and DenoiseReadCounts now run on tumor samples only.** Normal BAMs are only used for CollectAllelicCounts (required by ModelSegments). Normal QC (denoising, segmentation) belongs to the companion pipeline `gatk-cnv-pon`.

### Added

- **External gene list** (`assets/cancer_genes.txt`) — ~107 cancer genes previously hardcoded in the R script. Users can supply their own list via `--genes_list`.
- **`--genes_list` parameter** — Points to the gene list file. Defaults to `${projectDir}/assets/cancer_genes.txt`.
- **`--somatic_min_contig_length` parameter** — Configurable minimum contig length for PlotDenoisedCopyRatios and PlotModeledSegments (default: 46709983, chr22 in hg38).
- **`conf/modules.config`** — Centralized `publishDir` configuration and GATK 4.6.2.0 container overrides.
- **`conda` directive in `ANNOTATE_CNV`** — Enables execution without Docker.
- **Annotation fallback** — If the R script fails, an empty output with correct headers is created and the pipeline continues.
- **`tumor_name` field in meta map** — Supports multiple tumors per patient.
- **Execution reports** — timeline, report, trace, and dag enabled by default.

### Changed

- **R script moved to `bin/`** — `annotate_cnvs_complete_fixed.R` → `bin/annotate_cnvs.R`. No longer passed as process input.
- **R script accepts `--genes_list`** — Reads gene symbols from external file instead of hardcoded list.
- **R script comments in English** — Translated for international visibility.
- **GATK version standardized to 4.6.2.0** — All local modules overridden via `conf/modules.config`.
- **Tumor-only processing in steps 1-2** — CollectReadCounts and DenoiseReadCounts no longer process normal samples. Normal denoised/segmented results were computed but discarded; this wasted ~40% of compute time for those steps.
- **Robust double-key joins** — Tumor files joined by `meta.id` (unique per tumor), normal paired via `combine(by: 0)` on `patient_id`. Supports multiple tumors per patient.
- **FASTA documented as "main chroms only"** — Prominent comments in workflow, config, and schema.
- **Plot modules receive `--minimum-contig-length`** via `ext.args` in `conf/modules.config`.

### Fixed

- **Duplicate output channels in `ANNOTATE_CNV`** — `annotated` and `summary` both matched `*.txt`. Now each has a distinct pattern.

### Removed

- **`params.annotate_script`** — R script bundled in `bin/`.
- **Hardcoded gene list** — Replaced by `assets/cancer_genes.txt`.
- **Normal processing in CollectReadCounts/DenoiseReadCounts** — Moved to `gatk-cnv-pon` pipeline.

---

## [1.0.0] - 2025-XX-XX

### Added

- Initial release of the GATK CNV Somatic pipeline.
- Nextflow DSL2 implementation based on nf-core template.
- Tumor-normal paired CNV analysis using GATK4 best practices.
- Docker images: `florpio/cnv-references:hg38-v1.0`, `florpio/cnv-annotate-r:1.0`.
