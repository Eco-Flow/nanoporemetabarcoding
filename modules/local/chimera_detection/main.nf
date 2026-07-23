process CHIMERA_DETECTION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c86cbb145786bf5c24ea7fb13448da5f7d5cd124fd4403c1da5bc8fc60c2588/data':
        'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(db)

    output:
    tuple val(meta), path("${meta.id}.nonchimeric.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}.chimeras.tsv")     , emit: report
    tuple val("${task.process}"), val("blastn"), eval("blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def min_pid  = params.chimera_min_identity ?: 90    // Minimum % identity for a half's hit to count
    def min_qcov = params.chimera_min_coverage ?: 80    // Minimum % coverage of the half by its alignment
    def top_n    = params.chimera_top_n ?: 5            // Number of top reference hits to consider per half
    def margin   = params.chimera_bitscore_margin ?: 10 // Keep hits within this % of a half's best bitscore
    def is_compressed = fasta.getExtension() == "gz"
    def fasta_name    = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    # 1. Split every sequence into a 5' half (__H5) and a 3' half (__H3)
    awk '
        /^>/ { if (id != "") emit(); id = substr(\$1, 2); seq = ""; next }
        { seq = seq \$0 }
        END { if (id != "") emit() }
        function emit(   L, h) {
            L = length(seq); h = int(L / 2)
            if (h < 1 || L - h < 1) return          # too short to split, skip
            print ">" id "__H5\\n" substr(seq, 1, h)
            print ">" id "__H3\\n" substr(seq, h + 1)
        }
    ' ${fasta_name} > halves.fasta

    export BLASTDB=${db}
    DB=`find -L ./ -name "*.nal" | sed 's/\\.nal\$//'`
    if [ -z "\$DB" ]; then
        DB=`find -L ./ -name "*.nin" | sed 's/\\.nin\$//'`
    fi
    echo Using \$DB

    # 2. BLAST each half against the reference database, keeping the top N hits
    blastn \\
        -num_threads ${task.cpus} \\
        -db \$DB \\
        -query halves.fasta \\
        -outfmt '6 qseqid sseqid pident qcovhsp bitscore' \\
        -max_target_seqs ${top_n} \\
        -out halves.blast.tsv

    # 3. Build the set of plausible reference subjects for each half (hits passing the
    #    identity/coverage cut-offs and within 'margin' % of that half's best bitscore).
    #    A sequence is only flagged as chimeric when the two halves have NO reference in
    #    common - i.e. neither half's hits can explain the other. This tolerates the
    #    common case where redundant database entries for the same organism differ.
    awk -v pid=${min_pid} -v qcov=${min_qcov} -v margin=${margin} '
        function collect(id, half, out,    key, i, n, thr) {
            key = id SUBSEP half
            n = cnt[key]
            if (n == 0) return 0
            thr = best[key] * (1 - margin / 100)
            for (i = 1; i <= n; i++) if (bs[key, i] >= thr) out[subj[key, i]] = 1
            return 1
        }
        function join(arr,   s, k) { s = ""; for (k in arr) s = s (s ? "," : "") k; return s }
        {
            split(\$1, a, "__")
            half = a[length(a)]                     # H5 or H3
            id   = substr(\$1, 1, length(\$1) - length(half) - 2)
            if (\$3 + 0 >= pid && \$4 + 0 >= qcov) {
                key = id SUBSEP half
                n = ++cnt[key]
                subj[key, n] = \$2; bs[key, n] = \$5 + 0
                if (bs[key, n] > best[key]) best[key] = bs[key, n]
                ids[id] = 1
            }
        }
        END {
            print "seq_id\\th5_subjects\\th3_subjects\\tshared_subjects\\tchimeric"
            for (id in ids) {
                delete o5; delete o3
                h5 = collect(id, "H5", o5); h3 = collect(id, "H3", o3)
                nshared = 0
                for (s in o5) if (s in o3) nshared++
                flag = (h5 && h3 && nshared == 0) ? "yes" : "no"
                shared = ""
                for (s in o5) if (s in o3) shared = shared (shared ? "," : "") s
                print id "\\t" (h5 ? join(o5) : "NA") "\\t" (h3 ? join(o3) : "NA") "\\t" (nshared ? shared : "NA") "\\t" flag
            }
        }
    ' halves.blast.tsv | sort > ${prefix}.chimeras.tsv

    # 4. Write the non-chimeric sequences (everything not flagged "yes")
    awk -F'\\t' '\$5 == "yes" { print \$1 }' ${prefix}.chimeras.tsv > chimera_ids.txt
    awk '
        NR == FNR { drop[\$1] = 1; next }
        /^>/ { keep = !(substr(\$1, 2) in drop) }
        keep { print }
    ' chimera_ids.txt ${fasta_name} > ${prefix}.nonchimeric.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nonchimeric.fasta
    touch ${prefix}.chimeras.tsv
    """
}
