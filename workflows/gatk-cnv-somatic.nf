/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// nf-core modules
include { GATK4_PREPROCESSINTERVALS } from '../modules/nf-core/gatk4/preprocessintervals/main'
include { GATK4_COLLECTREADCOUNTS   } from '../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DENOISEREADCOUNTS   } from '../modules/nf-core/gatk4/denoisereadcounts/main'

// Local modules
include { COLLECTALLELICCOUNTS    } from '../modules/local/collectalleliccounts'
include { MODELSEGMENTS_SOMATIC   } from '../modules/local/modelsegments_somatic'
include { CALLCOPYRATIOSEGMENTS   } from '../modules/local/callcopyratiosegments'
include { PLOTDENOISEDCOPYRATIOS  } from '../modules/local/plotdenoisedcopyratios'
include { PLOTMODELEDSEGMENTS     } from '../modules/local/plotmodeledsegments'
include { ANNOTATE_CNV            } from '../modules/local/annotate_cnv'

// nf-core subworkflows
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GATK_CNV_SOMATIC {

    take:
    ch_samplesheet  // channel: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai,
                    //            tumor_counts, tumor_allelic, normal_allelic ]
                    // Any BAM/precomputed file may be [] (sentinel for "not provided"),
                    // validated upstream in PIPELINE_INITIALISATION.

    main:

    ch_versions = Channel.empty()

    // =========================================
    // Prepare reference channels
    // NOTE: FASTA must be main chromosomes only
    //       (hg38 without alt/decoy/HLA contigs)
    //       to match the Panel of Normals and
    //       preprocessed interval list.
    // =========================================

    ch_fasta       = Channel.value([ [:], file(params.fasta) ])
    ch_fai         = Channel.value([ [:], file(params.fai) ])
    ch_dict        = Channel.value([ [:], file(params.dict) ])
    ch_pon         = Channel.value([ [:], file(params.pon) ])
    ch_common_snps = file(params.common_snps)
    ch_common_snps_tbi = file(params.common_snps + ".tbi")

    // =========================================
    // 0. PREPARE INTERVALS
    //
    // Auto-detect input type:
    //   .interval_list → use directly (skip preprocessing)
    //   .bed           → run PreprocessIntervals
    //
    // IMPORTANT: The interval list must match the one
    // used to create the Panel of Normals. If you provide
    // a .bed file, ensure the same .bed + padding + bin_length
    // were used for PON creation.
    // =========================================

    def intervals_file = file(params.intervals)

    if (intervals_file.name.endsWith('.interval_list')) {

        log.info "Intervals: using pre-computed interval list (skipping PreprocessIntervals)"
        ch_intervals = intervals_file

    } else {

        log.info "Intervals: running PreprocessIntervals on ${intervals_file.name}"

        ch_intervals_input = params.intervals
            ? Channel.value([ [:], file(params.intervals) ])
            : Channel.value([ [:], [] ])
        ch_exclude_intervals = params.exclude_intervals
            ? Channel.value([ [:], file(params.exclude_intervals) ])
            : Channel.value([ [:], [] ])

        GATK4_PREPROCESSINTERVALS(
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_intervals_input,
            ch_exclude_intervals
        )
        ch_versions = ch_versions.mix(GATK4_PREPROCESSINTERVALS.out.versions)

        ch_intervals = GATK4_PREPROCESSINTERVALS.out.interval_list
            .map { meta, interval_list -> interval_list }
    }

    // =========================================
    // Prepare per-sample channels
    //
    // For each samplesheet row we derive three "task" channels — one per
    // process that consumes BAMs. Rows that already provide the precomputed
    // output (tumor_counts / tumor_allelic / normal_allelic) bypass the
    // corresponding process; rows that only have BAMs run it as before.
    // Mixed samplesheets are supported.
    //
    // Normal QC (denoising, segmentation) is handled by the gatk-cnv-pon
    // pipeline, not here.
    // =========================================

    // Per-row tumor view: [ tumor_meta, t_bam, t_bai, t_counts, t_allelic ]
    ch_tumor_row = ch_samplesheet
        .map { meta, t_bam, t_bai, n_bam, n_bai, t_counts, t_allelic, n_allelic ->
            def tumor_meta = meta + [
                sample_type: 'tumor',
                id:          "${meta.id}_tumor",
                patient_id:  meta.id,
                tumor_name:  "${meta.id}_tumor"
            ]
            [ tumor_meta, t_bam, t_bai, t_counts, t_allelic ]
        }

    // Per-row normal view: [ normal_meta, n_bam, n_bai, n_allelic ]
    ch_normal_row = ch_samplesheet
        .map { meta, t_bam, t_bai, n_bam, n_bai, t_counts, t_allelic, n_allelic ->
            def normal_meta = meta + [
                sample_type: 'normal',
                id:          "${meta.id}_normal",
                patient_id:  meta.id
            ]
            [ normal_meta, n_bam, n_bai, n_allelic ]
        }

    // =========================================
    // 1. COLLECT READ COUNTS (tumor only)
    //    Skip rows where tumor_counts is already provided.
    // =========================================

    ch_tumor_row
        .branch { meta, bam, bai, counts, allelic ->
            precomputed: counts as boolean
            from_bam:    !counts
        }
        .set { ch_tumor_counts_branch }

    ch_tumor_for_collect = ch_tumor_counts_branch.from_bam
        .map { meta, bam, bai, counts, allelic -> [ meta, bam, bai, ch_intervals ] }

    GATK4_COLLECTREADCOUNTS(
        ch_tumor_for_collect,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions.first())

    // Unified tumor read counts channel: [ meta, hdf5 ]
    ch_tumor_counts = GATK4_COLLECTREADCOUNTS.out.hdf5
        .mix(
            ch_tumor_counts_branch.precomputed
                .map { meta, bam, bai, counts, allelic -> [ meta, counts ] }
        )

    // =========================================
    // 2. DENOISE READ COUNTS (tumor only)
    //    Runs for every tumor (precomputed or freshly collected).
    // =========================================

    GATK4_DENOISEREADCOUNTS(
        ch_tumor_counts,
        ch_pon
    )
    ch_versions = ch_versions.mix(GATK4_DENOISEREADCOUNTS.out.versions.first())

    // =========================================
    // 3. COLLECT ALLELIC COUNTS (tumor + normal)
    //    Skip rows where the corresponding allelic_counts file is already
    //    provided. Tumor and normal are branched independently so a row can
    //    have e.g. tumor_allelic precomputed but normal_bam fresh.
    //
    //    Both outputs are needed by ModelSegments:
    //      tumor allelic  → --allelic-counts
    //      normal allelic → --normal-allelic-counts
    // =========================================

    ch_tumor_row
        .branch { meta, bam, bai, counts, allelic ->
            precomputed: allelic as boolean
            from_bam:    !allelic
        }
        .set { ch_tumor_allelic_branch }

    ch_normal_row
        .branch { meta, bam, bai, allelic ->
            precomputed: allelic as boolean
            from_bam:    !allelic
        }
        .set { ch_normal_allelic_branch }

    ch_bams_for_allelic = ch_tumor_allelic_branch.from_bam
        .map { meta, bam, bai, counts, allelic -> [ meta, bam, bai ] }
        .mix(
            ch_normal_allelic_branch.from_bam
                .map { meta, bam, bai, allelic -> [ meta, bam, bai ] }
        )

    COLLECTALLELICCOUNTS(
        ch_bams_for_allelic,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_common_snps,
        ch_common_snps_tbi
    )
    ch_versions = ch_versions.mix(COLLECTALLELICCOUNTS.out.versions.first())

    // Unified allelic counts channel: [ meta, allelic_tsv ]
    // (tumor + normal mixed, downstream filtered by meta.sample_type)
    ch_all_allelic = COLLECTALLELICCOUNTS.out.allelic_counts
        .mix(
            ch_tumor_allelic_branch.precomputed
                .map { meta, bam, bai, counts, allelic -> [ meta, allelic ] }
        )
        .mix(
            ch_normal_allelic_branch.precomputed
                .map { meta, bam, bai, allelic -> [ meta, allelic ] }
        )

    // =========================================
    // 4. MODEL SEGMENTS (tumor with matched normal)
    //
    // Join strategy:
    //   - Tumor files joined by meta.id (unique per tumor,
    //     supports multiple tumors per patient)
    //   - Normal paired via combine(by: 0) on patient_id
    //     (1 normal broadcast to N tumors)
    // =========================================

    ch_tumor_denoised = GATK4_DENOISEREADCOUNTS.out.denoised
        .map { meta, file -> [ meta.id, meta, file ] }

    ch_tumor_allelic = ch_all_allelic
        .filter { meta, file -> meta.sample_type == 'tumor' }
        .map { meta, file -> [ meta.id, file ] }

    ch_normal_allelic = ch_all_allelic
        .filter { meta, file -> meta.sample_type == 'normal' }
        .map { meta, file -> [ meta.patient_id, file ] }

    // Step 1: join tumor denoised + tumor allelic by meta.id (1:1)
    // Step 2: re-key by patient_id, combine with normal (1:N)
    ch_model_segments_input = ch_tumor_denoised
        .join(ch_tumor_allelic)
        .map { tumor_id, meta, denoised, tumor_allelic ->
            [ meta.patient_id, meta, denoised, tumor_allelic ]
        }
        .combine(ch_normal_allelic, by: 0)
        .map { patient_id, meta, denoised, tumor_allelic, normal_allelic ->
            [ meta, denoised, tumor_allelic, normal_allelic ]
        }

    MODELSEGMENTS_SOMATIC(
        ch_model_segments_input
    )
    ch_versions = ch_versions.mix(MODELSEGMENTS_SOMATIC.out.versions.first())

    // =========================================
    // 5. CALL COPY RATIO SEGMENTS
    // =========================================

    CALLCOPYRATIOSEGMENTS(
        MODELSEGMENTS_SOMATIC.out.cr_segments
    )
    ch_versions = ch_versions.mix(CALLCOPYRATIOSEGMENTS.out.versions.first())

    // =========================================
    // 6. PLOTS (optional)
    //    --minimum-contig-length is passed via
    //    ext.args in conf/modules.config using
    //    params.somatic_min_contig_length
    // =========================================

    if (!params.skip_plots) {

        ch_tumor_standardized = GATK4_DENOISEREADCOUNTS.out.standardized
            .map { meta, file -> [ meta.id, file ] }

        ch_plot_denoised_input = ch_tumor_denoised
            .join(ch_tumor_standardized)
            .map { tumor_id, meta, denoised, standardized ->
                [ meta, standardized, denoised ]
            }

        PLOTDENOISEDCOPYRATIOS(
            ch_plot_denoised_input,
            ch_dict
        )
        ch_versions = ch_versions.mix(PLOTDENOISEDCOPYRATIOS.out.versions.first())

        ch_plot_modeled_input = ch_tumor_denoised
            .join(
                MODELSEGMENTS_SOMATIC.out.hets
                    .map { meta, file -> [ meta.id, file ] }
            )
            .join(
                MODELSEGMENTS_SOMATIC.out.segments
                    .map { meta, file -> [ meta.id, file ] }
            )
            .map { tumor_id, meta, denoised, hets, segments ->
                [ meta, denoised, hets, segments ]
            }

        PLOTMODELEDSEGMENTS(
            ch_plot_modeled_input,
            ch_dict
        )
        ch_versions = ch_versions.mix(PLOTMODELEDSEGMENTS.out.versions.first())
    }

    // =========================================
    // 7. ANNOTATE CNVs (optional)
    //    R script is in bin/ (auto-available in $PATH)
    //    Has fallback: creates empty output if R fails
    // =========================================

    if (!params.skip_annotation) {

        ch_annotate_input = MODELSEGMENTS_SOMATIC.out.segments
            .join(MODELSEGMENTS_SOMATIC.out.cr_segments, by: 0)
            .map { meta, segments, cr_seg ->
                [ meta, segments ]
            }

        ANNOTATE_CNV(
            ch_annotate_input,
            file(params.mane_annotation),
            file(params.genes_list)
        )
        ch_versions = ch_versions.mix(ANNOTATE_CNV.out.versions.first())
    }

    // =========================================
    // Collate and save software versions
    // =========================================
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'gatk_cnv_somatic_software_versions.yml',
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

    emit:
    segments         = MODELSEGMENTS_SOMATIC.out.segments
    called_segments  = CALLCOPYRATIOSEGMENTS.out.called_segments
    versions         = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
