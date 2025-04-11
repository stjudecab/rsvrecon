process IGVTOOLS_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::igvtools=2.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/igvtools:2.14.1--hdfd78af_0' :
        'biocontainers/igvtools:2.14.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(fasta)

    output:
    tuple val(meta), path("*.wig"), emit: count
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export LC_ALL=C
    export _JAVA_OPTIONS="-Duser.home=\$(pwd)"

    igvtools count \\
        $args \\
        --bases $bam \\
        "${prefix}.coverage.wig" \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igvtools: \$(igvtools version | grep -o "Version [0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed 's/Version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.coverage.wig"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igvtools: \$(igvtools version | grep -o "Version [0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed 's/Version //')
    END_VERSIONS
    """
}
