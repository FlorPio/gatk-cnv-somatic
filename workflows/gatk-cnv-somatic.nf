/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// nf-core modules
include { GATK4_COLLECTREADCOUNTS } from '../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DENOISEREADCOUNTS } from '../modules/nf-core/gatk4/denoisereadcounts/main'

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
    ch_samplesheet  // channel: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]

    main:

    ch_versions = Channel.empty()

    // =========================================
    // Prepare reference channels
    // =========================================
    
    ch_fasta       = Channel.value([ [:], file(params.fasta) ])
    ch_fai         = Channel.value([ [:], file(params.fai) ])
    ch_dict        = Channel.value([ [:], file(params.dict) ])
    ch_pon         = Channel.value([ [:], file(params.pon) ])
    ch_intervals   = file(params.intervals)
    ch_common_snps = file(params.common_snps)
    ch_common_snps_tbi = file(params.common_snps + ".tbi")

    // =========================================
    // Prepare BAM channels (tumor + normal combined)
    // =========================================
    
    ch_all_bams = ch_samplesheet
        .flatMap { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            def tumor_meta  = meta + [sample_type: 'tumor',  id: "${meta.id}_tumor", patient_id: meta.id]
            def normal_meta = meta + [sample_type: 'normal', id: "${meta.id}_normal", patient_id: meta.id]
            [
                [ tumor_meta,  tumor_bam,  tumor_bai ],
                [ normal_meta, normal_bam, normal_bai ]
            ]
        }

    ch_bams_with_intervals = ch_all_bams
        .map { meta, bam, bai -> [ meta, bam, bai, ch_intervals ] }

    // =========================================
    // 1. COLLECT READ COUNTS (tumor + normal)
    // =========================================
    
    GATK4_COLLECTREADCOUNTS(
        ch_bams_with_intervals,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions.first())

    // =========================================
    // 2. DENOISE READ COUNTS (tumor + normal)
    // =========================================
    
    GATK4_DENOISEREADCOUNTS(
        GATK4_COLLECTREADCOUNTS.out.hdf5,
        ch_pon
    )
    ch_versions = ch_versions.mix(GATK4_DENOISEREADCOUNTS.out.versions.first())

    // =========================================
    // 3. COLLECT ALLELIC COUNTS (tumor + normal)
    // =========================================
    
    COLLECTALLELICCOUNTS(
        ch_all_bams,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_common_snps,
        ch_common_snps_tbi
    )
    ch_versions = ch_versions.mix(COLLECTALLELICCOUNTS.out.versions.first())

    // =========================================
    // 4. MODEL SEGMENTS (tumor with matched normal)
    // =========================================
    
    ch_tumor_denoised = GATK4_DENOISEREADCOUNTS.out.denoised
        .filter { meta, file -> meta.sample_type == 'tumor' }
        .map { meta, file -> [ meta.patient_id, meta, file ] }

    ch_tumor_standardized = GATK4_DENOISEREADCOUNTS.out.standardized
        .filter { meta, file -> meta.sample_type == 'tumor' }
        .map { meta, file -> [ meta.patient_id, file ] }

    ch_tumor_allelic = COLLECTALLELICCOUNTS.out.allelic_counts
        .filter { meta, file -> meta.sample_type == 'tumor' }
        .map { meta, file -> [ meta.patient_id, file ] }

    ch_normal_allelic = COLLECTALLELICCOUNTS.out.allelic_counts
        .filter { meta, file -> meta.sample_type == 'normal' }
        .map { meta, file -> [ meta.patient_id, file ] }

    ch_model_segments_input = ch_tumor_denoised
        .join(ch_tumor_allelic)
        .join(ch_normal_allelic)
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
    // =========================================
    
    if (!params.skip_plots) {
        
        ch_plot_denoised_input = ch_tumor_denoised
            .join(ch_tumor_standardized)
            .map { patient_id, meta, denoised, standardized ->
                [ meta, standardized, denoised ]
            }

        PLOTDENOISEDCOPYRATIOS(
            ch_plot_denoised_input,
            ch_dict
        )
        ch_versions = ch_versions.mix(PLOTDENOISEDCOPYRATIOS.out.versions.first())

        ch_plot_modeled_input = GATK4_DENOISEREADCOUNTS.out.denoised
            .filter { meta, file -> meta.sample_type == 'tumor' }
            .map { meta, file -> [ meta.patient_id, meta, file ] }
            .join(
                MODELSEGMENTS_SOMATIC.out.hets
                    .map { meta, file -> [ meta.patient_id, file ] }
            )
            .join(
                MODELSEGMENTS_SOMATIC.out.segments
                    .map { meta, file -> [ meta.patient_id, file ] }
            )
            .map { patient_id, meta, denoised, hets, segments ->
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
            file(params.annotate_script)
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
