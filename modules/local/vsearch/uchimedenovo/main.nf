process VSEARCH_UCHIMEDENOVO {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1' :
        'biocontainers/vsearch:2.21.1--hf1761c0_1' }"

    input:
    tuple val(meta), path(fasta), path(counts)

    output:
    tuple val(meta), path("${meta.id}.nonchimeras.fasta"), emit: nonchimeras
    tuple val(meta), path("${meta.id}.chimeras.fasta")   , emit: chimeras
    tuple val(meta), path("${meta.id}.uchime.txt")       , emit: report   , optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '--uchime3_denovo'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz"
    def fasta_name    = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    # De novo chimera detection relies on abundance skew (chimeras are rarer than
    # their parents), which vsearch reads from a ';size=N' annotation on each header.
    # Annotate abundances from the per-ASV read counts (CSV: 'seq_id,count'); any
    # sequence without a matching count falls back to size=1.
    if [ -s "${counts}" ]; then
        awk '
            BEGIN { FS = "," }
            FNR == NR {
                id = \$1; c = \$2
                gsub(/^[ \\t]+|[ \\t]+\$/, "", id)
                gsub(/^[ \\t]+|[ \\t]+\$/, "", c)
                if (id != "") size[id] = c
                next
            }
            /^>/ {
                split(substr(\$0, 2), p, /[ \\t]/)
                name = p[1]
                s = (name in size) ? size[name] : 1
                print ">" name ";size=" s
                next
            }
            { print }
        ' ${counts} ${fasta_name} > sized.fasta
    else
        awk '/^>/ && \$0 !~ /;size=/ { print \$1 ";size=1"; next } { print }' ${fasta_name} > sized.fasta
    fi

    vsearch \\
        ${args} sized.fasta \\
        --sizein \\
        --nonchimeras ${prefix}.nonchimeras.fasta \\
        --chimeras ${prefix}.chimeras.fasta \\
        --uchimeout ${prefix}.uchime.txt \\
        --threads ${task.cpus}

    # Guarantee the emitted files exist even when vsearch finds nothing in a class
    touch ${prefix}.nonchimeras.fasta ${prefix}.chimeras.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nonchimeras.fasta
    touch ${prefix}.chimeras.fasta
    touch ${prefix}.uchime.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
