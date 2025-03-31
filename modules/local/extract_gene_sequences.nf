process EXTRACT_GENE_SEQUENCES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9.12 bioconda::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' :
        'quay.io/biocontainers/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' }"

    input:
    tuple val(meta), path(assembled_genome), path(reference_genome), path(gff_file)
    val cds_names

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in stjudecab/rsvrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cds_list = cds_names instanceof List ? cds_names.join(',') : cds_names
    """
    extract_gene_sequences.py \\
        -a $assembled_genome \\
        -r $reference_genome \\
        -g $gff_file \\
        -c $cds_list \\
        -o "${prefix}.ggene.fasta" \\
        --log "${prefix}.log" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
