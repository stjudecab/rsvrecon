process PARSE_GGENE_BLASTN {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9.12 conda-forge::pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' :
        'quay.io/biocontainers/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' }"

    input:
    tuple val(meta), path(blastn_out)

    output:
    tuple val(meta), path("*.ggene_genotype.txt"), emit: genotype
    tuple val(meta), path("*.log")               , emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in stjudecab/rsvrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_ggene_blastn.py \\
        --blast_out ${blastn_out} \\
        --output ${prefix}.ggene_genotype.txt \\
        --seq_id ${meta.id} \\
        --log "${prefix}.log" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
