process AMPLICON_SORTER {
    tag "$meta.id"
    label 'process_medium'

    container 'community.wave.seqera.io/library/pip_biopython_edlib_matplotlib:2788dbec4d078d0a'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}/consensusfile.fasta")   , emit: consensus, optional: true
    tuple val(meta), path("${meta.id}/results.txt")           , emit: results_txt, optional: true
    tuple val(meta), path("${meta.id}/results.csv")           , emit: results_csv, optional: true
    tuple val(meta), path("${meta.id}/*_[0-9]*_[0-9]*.fasta") , emit: fastas, optional: true // Numbers only as there are unisgned sequences with _nongroup_ and _unique_ labels
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 ${projectDir}/bin/amplicon_sorter.py \\
        -i $reads \\
        -o $prefix \\
        -np $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(python3 --version)
    END_VERSIONS
    """

}
