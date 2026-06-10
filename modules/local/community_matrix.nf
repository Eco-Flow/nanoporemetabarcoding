process COMMUNITY_MATRIX {
    tag "$meta.id"
    label 'process_single'

    //conda
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-argparse_r-dplyr_r-tidyr_r-vegan:cbe9ef99325d05b4':
        'community.wave.seqera.io/library/r-argparse_r-dplyr_r-tidyr_r-vegan:d6b6b06cd1b14ebb' }"

    input:
    tuple val(meta), path(table)

    output:
    tuple val(meta), path("*.csv")                                                                                                      , emit: table

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    Rscript ${projectDir}/bin/community_matrix.R \\
        $table \\
        --prefix ${prefix} \\
        $args
    """
}