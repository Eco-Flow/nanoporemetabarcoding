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
    path(sql_db)

    //output:
    //tuple val(meta), path("${meta.id}_best_hit.txt")   , emit: best_hit
    //path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript ${projectDir}/bin/assign_taxonomy.R $blast_hits \\
    --sql_db $sql_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version)
    END_VERSIONS
    """
}
