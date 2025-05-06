process ASSEMBLE_SEQUENCE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9.12 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' :
        'quay.io/biocontainers/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' }"

    input:
    tuple val(meta), path(count), path(fasta)
    val igv_cutoff

    output:
    tuple val(meta), path("*.consensus.fasta"), emit: consensus
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in stjudecab/rsvrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    assemble_sequence.py \\
        --wig $count \\
        --ref-fasta $fasta \\
        --igv-cutoff $igv_cutoff \\
        --out "${prefix}.consensus.fasta" \\
        --id "$meta.id" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

