process COLLECTALLELICCOUNTS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::gatk4=4.4.0.0"
    container "docker.io/broadinstitute/gatk:4.4.0.0"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path(common_snps)
    path(common_snps_tbi)

    output:
    tuple val(meta), path("*.allelicCounts.tsv"), emit: allelic_counts
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK CollectAllelicCounts] Available memory not known - defaulting to 3GB.')
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CollectAllelicCounts \\
        -I ${bam} \\
        -L ${common_snps} \\
        -R ${fasta} \\
        -O ${prefix}.allelicCounts.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.allelicCounts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
