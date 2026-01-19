#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FlorPio/gatk-cnv-somatic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/FlorPio/gatk-cnv-somatic
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK_CNV_SOMATIC  } from './workflows/gatk-cnv-somatic'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_gatk-cnv-somatic_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_gatk-cnv-somatic_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow FLORPIO_GATK_CNV_SOMATIC {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    GATK_CNV_SOMATIC (
        samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    FLORPIO_GATK_CNV_SOMATIC (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
