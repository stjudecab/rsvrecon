process KMA_MAP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kma=1.4.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.15--h577a1d6_1' :
        'biocontainers/kma:1.4.15--h577a1d6_1' }"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("${prefix}.res"),      emit: res      // Results overview
    tuple val(meta), path("${prefix}.fsa"),      emit: fsa      // Consensus sequences
    tuple val(meta), path("${prefix}.aln"),      emit: aln      // Consensus alignments
    tuple val(meta), path("${prefix}.frag.gz"),  emit: frag     // Read mapping information
    tuple val(meta), path("${prefix}.mat.gz"),   optional:true, emit: matrix   // Base counts (only if -matrix is enabled)
    tuple val(meta), path("${prefix}.log"),      emit: log      // Log file
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def args       = task.ext.args ?: ''
    def index_name = meta2.id ?: "index"  // if index has no id, then use 'index' as default. So, ensure meta2 contains the name of index

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
    kma ${read_command} \\
        -o ${prefix} \\
        -t_db ${index}/${index_name} \\
        $args \\
        2> ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma -v 2>&1) | sed 's/^KMA-//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${fasta.baseName}"
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
