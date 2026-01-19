process PLOTDENOISEDCOPYRATIOS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::gatk4=4.4.0.0"
    container "docker.io/broadinstitute/gatk:4.4.0.0"

    input:
    tuple val(meta), path(standardized_cr), path(denoised_cr)
    tuple val(meta2), path(dict)

    output:
    tuple val(meta), path("*.denoised.png"),      emit: denoised_plot
    tuple val(meta), path("*_standardized.png"),  emit: standardized_plot, optional: true
    tuple val(meta), path("*.denoised.png"),      emit: plot_png, optional: true
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK PlotDenoisedCopyRatios] Available memory not known - defaulting to 3GB.')
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        PlotDenoisedCopyRatios \\
        --standardized-copy-ratios ${standardized_cr} \\
        --denoised-copy-ratios ${denoised_cr} \\
        --sequence-dictionary ${dict} \\
        --output . \\
        --output-prefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_denoised.png
    touch ${prefix}_standardized.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
