process MODELSEGMENTS_SOMATIC {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::gatk4=4.4.0.0"
    container "docker.io/broadinstitute/gatk:4.6.2.0"

    input:
    tuple val(meta), path(denoised_cr), path(tumor_allelic), path(normal_allelic)

    output:
    tuple val(meta), path("*.modelFinal.seg"),  emit: segments
    tuple val(meta), path("*.cr.seg"),          emit: cr_segments
    tuple val(meta), path("*.hets.tsv"),        emit: hets
    tuple val(meta), path("*.modelBegin.seg"),  emit: model_begin, optional: true
    tuple val(meta), path("*.modelFinal.af.param"), emit: af_param, optional: true
    tuple val(meta), path("*.modelFinal.cr.param"), emit: cr_param, optional: true
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK ModelSegments] Available memory not known - defaulting to 3GB.')
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ModelSegments \\
        --denoised-copy-ratios ${denoised_cr} \\
        --allelic-counts ${tumor_allelic} \\
        --normal-allelic-counts ${normal_allelic} \\
        --output-prefix ${prefix} \\
        -O . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.modelFinal.seg
    touch ${prefix}.cr.seg
    touch ${prefix}.hets.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
