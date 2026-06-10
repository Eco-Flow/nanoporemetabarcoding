process PLOT_TAXONOMY {
    tag "$meta.id"
    label 'process_single'

    //conda
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-argparse_r-dplyr_r-ggplot2_r-tidyr:ee5d96cb9bb55dfa' :
        'community.wave.seqera.io/library/r-argparse_r-dplyr_r-ggplot2_r-tidyr:b7b1d28a987c6231' }"

    input:
    tuple val(meta), path(tables, stageAs: 'table_????')

    output:
    tuple val(meta), path("*.png")                                                                                                      , emit: plot
    tuple val("${task.process}"), val("R"), eval('Rscript -e "cat(as.character(getRversion()))"')                                       , topic: versions, emit: versions_r
    tuple val("${task.process}"), val("ggplot2"), eval('Rscript -e "cat(as.character(packageVersion(\'ggplot2\')))"')                   , topic: versions, emit: versions_ggplot2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sed '2,\${/sample_name/d;}' $tables > ${prefix}_join_table.csv

    Rscript ${projectDir}/bin/plot_asv_per_phylum.R \\
        ${prefix}_join_table.csv \\
        $args
    """
}
