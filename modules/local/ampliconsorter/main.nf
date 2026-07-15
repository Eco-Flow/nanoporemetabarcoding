process AMPLICON_SORTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/pip_biopython_edlib_matplotlib:2788dbec4d078d0a'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}/consensusfile.fasta")   , emit: consensus, optional: true
    tuple val(meta), path("${meta.id}/results.txt")           , emit: results_txt, optional: true
    tuple val(meta), path("${meta.id}/results.csv")           , emit: results_csv, optional: true
    tuple val(meta), path("${meta.id}/*/*_[0-9]*_[0-9]*.fasta") , emit: fastas, optional: true // Numbers only as there are unisgned sequences with _nongroup_ and _unique_ labels
    tuple val("${task.process}"), val("python"), eval("python3 --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val("amplicon_sorter"), eval("python3 ${projectDir}/bin/amplicon_sorter.py --version"), topic: versions, emit: versions_amplicon_sorter
    tuple val("${task.process}"), val("edlib"), eval("pip freeze | grep edlib | cut -d'=' -f3"), topic: versions, emit: versions_edlib
    tuple val("${task.process}"), val("biopython"), eval("pip freeze | grep biopython | cut -d'=' -f3"), topic: versions, emit: versions_biopython
    tuple val("${task.process}"), val("matplotlib"), eval("pip freeze | grep matplotlib | cut -d'=' -f3"), topic: versions, emit: versions_matplotlib

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

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    mkdir -p ${prefix}
    touch ${prefix}/consensusfile.fasta
    touch ${prefix}/results.txt
    touch ${prefix}/results.csv
    """
}
