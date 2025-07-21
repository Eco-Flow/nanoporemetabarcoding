process BEST_HIT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=8.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    input:
    tuple val(meta), path(blast_hits)

    output:
    tuple val(meta), path("${meta.id}_best_hit.txt")   , emit: best_hit
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sort -k1,1 -k12,12gr -k11,11g -k3,3gr $blast_hits | awk '!seen[\$1]++' > ${prefix}_best_hit.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version)
    END_VERSIONS
    """
}
