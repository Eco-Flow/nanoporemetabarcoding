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
    tuple val("${task.process}"), val("R"), eval("R --version | head -n1 | cut -d' ' -f3"), topic: versions, emit: versions_r
    tuple val("${task.process}"), val("taxonomizr"), eval("Rscript -e 'packageVersion(\"taxonomizr\")' | sed 's/.*\\s//; s/[\\[\\]']//g'"), topic: versions, emit: versions_taxonomizr

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
