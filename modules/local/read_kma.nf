process READ_KMA {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9.12 conda-forge::biopython=1.79 conda-forge::pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' :
        'quay.io/biocontainers/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' }"

    input:
    tuple val(meta), path(kma_res)
    path fasta
    path gff
    path meta_file

    output:
    tuple val(meta), path("*.fasta"), env(REF_SUBTYPE), emit: fasta
    tuple val(meta), path("*.gff"),   env(REF_SUBTYPE), emit: gff
    tuple val(meta), path("*.vars")                   , emit: annotation
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in stjudecab/rsvrecon/bin/
    def args = task.ext.args ?: ''
    """
    read_kma.py \\
        --input $kma_res \\
        --ref-fasta $fasta \\
        --ref-gff $gff \\
        --ref-metadata $meta_file \\
        --output . \\
        $args \\
        > env.vars

    source env.vars

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

