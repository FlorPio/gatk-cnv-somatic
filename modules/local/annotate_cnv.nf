process ANNOTATE_CNV {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::r-base=4.3 conda-forge::r-optparse conda-forge::r-dplyr conda-forge::r-tidyr bioconda::bioconductor-genomicranges"
    container "docker.io/florpio/cnv-annotate-r:1.0"

    input:
    tuple val(meta), path(model_segments)
    path(mane_annotation)
    path(genes_list)

    output:
    tuple val(meta), path("*_annotated.txt"), emit: annotated
    tuple val(meta), path("*_purity.txt"),    emit: purity, optional: true
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotate_cnvs.R \\
        --input_dir . \\
        --mane ${mane_annotation} \\
        --genes_list ${genes_list} \\
        --output_dir . \\
        ${args}

    # Rename outputs to include sample prefix
    if [ -f "CNVs_annotated_all_samples.txt" ]; then
        mv CNVs_annotated_all_samples.txt ${prefix}_annotated.txt
    fi
    if [ -f "tumor_purity_estimates.txt" ]; then
        mv tumor_purity_estimates.txt ${prefix}_purity.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}

