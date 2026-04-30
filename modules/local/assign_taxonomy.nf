process ASSIGN_TAXONOMY {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::r-argparse=1.0.1 bioconda::r-taxonomizr=0.7.1 conda-forge::r-tidyverse=2.0.0"
    container 'quay.io/fduarte001/assign_taxa:v1.0'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/r-argparse:1.0.1--py36r3.3.2_0' :
    //    'quay.io/biocontainers/r-argparse:1.0.1--py36r3.3.2_0' }"

    input:
    tuple val(meta), path(blast_hits)
    tuple val(meta), path(read_counts)
    path(sql_db)

    output:
    tuple val(meta), path("*.csv")   , emit: tax_csv
    tuple val("${task.process}"), val("R"), eval('Rscript -e "cat(as.character(getRversion()))"'), topic: versions, emit: versions_r
    tuple val("${task.process}"), val("taxonomizr"), eval('Rscript -e "cat(as.character(packageVersion(\'taxonomizr\')))"'), topic: versions, emit: versions_taxonomizr
    tuple val("${task.process}"), val("argparse"), eval('Rscript -e "cat(as.character(packageVersion(\'argparse\')))"'), topic: versions, emit: versions_argparse
    tuple val("${task.process}"), val("dplyr"), eval('Rscript -e "cat(as.character(packageVersion(\'dplyr\')))"'), topic: versions, emit: versions_dplyr
    tuple val("${task.process}"), val("tidyr"), eval('Rscript -e "cat(as.character(packageVersion(\'tidyr\')))"'), topic: versions, emit: versions_tidyr
    tuple val("${task.process}"), val("stringr"), eval('Rscript -e "cat(as.character(packageVersion(\'stringr\')))"'), topic: versions, emit: versions_stringr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript ${projectDir}/bin/assign_taxonomy.R \\
    $blast_hits \\
    $read_counts \\
    --sql_db $sql_db \\
    $args

    """
}
