process KMA_MAP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kma=1.4.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.15--h577a1d6_1' :
        'biocontainers/kma:1.4.15--h577a1d6_1' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.res"),      emit: res      // Results overview
    tuple val(meta), path("*.fsa"),      emit: fsa      // Consensus sequences
    tuple val(meta), path("*.aln"),      emit: aln      // Consensus alignments
    tuple val(meta), path("*.frag.gz"),  emit: frag     // Read mapping information
    tuple val(meta), path("*.mat.gz"),   optional:true, emit: matrix   // Base counts (only if -matrix is enabled)
    tuple val(meta), path("*.log"),      emit: log      // Log file
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "$meta.id"

    // Handle different read formats
    def read_command = ""
    if ( task.ext.args?.contains('-int') ) {
        // Interleaved
        read_command = "-int ${reads}"
    } else if (meta.single_end) {
        // Single-end
        read_command = "-i ${reads}"
    } else {
        // Pair-end
        read_command = "-ipe ${reads[0]} ${reads[1]}"
    }

    """
    # Determine index base name by looking for required index files
    INDEX_BASE=""
    INDEX_FILES=\$(ls ${index}/*.seq.b 2>/dev/null || ls ${index}/*/*.seq.b 2>/dev/null || echo "")

    if [ -n "\$INDEX_FILES" ]; then
        # Extract the base name from the first matching index file
        INDEX_FILE=\$(echo "\$INDEX_FILES" | head -n 1)
        INDEX_BASE=\${INDEX_FILE%.seq.b}
        echo "Using index base: \$INDEX_BASE"
    else
        # If no *.seq.b files found, try to check if the index itself is the base name
        if [ -f "${index}.seq.b" ]; then
            INDEX_BASE="${index}"
            echo "Using index base: \$INDEX_BASE"
        else
            echo "Error: Could not find proper KMA index files" >&2
            exit 1
        fi
    fi

    # KMA will return nonezero code even no errors, so we add `|| true`
    kma ${read_command} \\
        -o ${prefix} \\
        -t_db \$INDEX_BASE \\
        $args \\
        2> ${prefix}.log || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma -v 2>&1) | sed 's/^KMA-//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.res
    touch ${prefix}.fsa
    touch ${prefix}.aln
    touch ${prefix}.frag.gz
    touch ${prefix}.log
    touch ${prefix}.mat.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma -v 2>&1) | sed 's/^KMA-//')
    END_VERSIONS
    """
}
