process ANNOTATE_CNV {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/florpio/cnv-annotate-r:1.0"

    input:
    tuple val(meta), path(model_segments)
    path(mane_annotation)
    path(rscript)

    output:
    tuple val(meta), path("*.txt"),     emit: annotated
    tuple val(meta), path("*.txt"),       emit: summary, optional: true
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript ${rscript} \\
        --input_dir . \\
        --mane ${mane_annotation} \\
        --output_dir . \\
        ${args}

    # Rename output to include sample prefix
    for f in *.txt; do
        if [ -f "\$f" ] && [[ "\$f" != "${prefix}"* ]]; then
            mv "\$f" "${prefix}_annotated.tsv" 2>/dev/null || true
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}
