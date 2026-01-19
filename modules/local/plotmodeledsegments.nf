process PLOTMODELEDSEGMENTS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::gatk4=4.4.0.0"
    container "docker.io/broadinstitute/gatk:4.4.0.0"

    input:
    tuple val(meta), path(denoised_cr), path(hets), path(segments)
    tuple val(meta2), path(dict)

    output:
    tuple val(meta), path("*.modeled.png"), emit: modeled_plot
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK PlotModeledSegments] Available memory not known - defaulting to 3GB.')
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        PlotModeledSegments \\
        --denoised-copy-ratios ${denoised_cr} \\
        --allelic-counts ${hets} \\
        --segments ${segments} \\
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
    touch ${prefix}.modeled.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
