process COMMUNITY_MATRIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-argparse_r-dplyr_r-tidyr_r-vegan:cbe9ef99325d05b4':
        'community.wave.seqera.io/library/r-argparse_r-dplyr_r-tidyr_r-vegan:d6b6b06cd1b14ebb' }"

    input:
    tuple val(meta), path(table)

    output:
    tuple val(meta), path("*.csv")
    tuple val("${task.process}"), val('argparse'), eval('Rscript -e "cat(as.character(packageVersion(\'argparse\')))"'), topic: versions, emit: versions_argparse
    tuple val("${task.process}"), val('dplyr'), eval('Rscript -e "cat(as.character(packageVersion(\'dplyr\')))"'), topic: versions, emit: versions_dplyr
    tuple val("${task.process}"), val('tidyr'), eval('Rscript -e "cat(as.character(packageVersion(\'tidyr\')))"'), topic: versions, emit: versions_tidyr
    tuple val("${task.process}"), val('vegan'), eval('Rscript -e "cat(as.character(packageVersion(\'vegan\')))"'), topic: versions, emit: versions_vegan

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

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.csv
    """
}
