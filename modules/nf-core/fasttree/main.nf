process FASTTREE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fasttree:2.1.10--h516909a_4' :
        'biocontainers/fasttree:2.1.10--h516909a_4' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.tre"), emit: phylogeny
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$meta.id"
    """
    fasttree \\
        $args \\
        -log fasttree_phylogeny.tre.log \\
        -nt $alignment \\
        > "${prefix}.fasttree_phylogeny.tre"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: \$(fasttree -help 2>&1 | head -1  | sed 's/^FastTree \\([0-9\\.]*\\) .*\$/\\1/')
    END_VERSIONS
    """
}
