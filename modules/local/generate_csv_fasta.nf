process GENERATE_CSV_FASTA {
    tag "GENERATE_CSV_FASTA"
    label 'process_low'

    conda "conda-forge::python=3.9.12 conda-forge::biopython=1.79 conda-forge::pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' :
        'quay.io/biocontainers/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' }"

    input:
    path(manifest, stageAs: "manifest/*")
    path(version , stageAs: "manifest/*")
    path(reference_dir)
    val igv_cutoff

    output:
    path "Report.csv", emit: report
    tuple path("tree_sequences_A.fasta", optional: true),
        path("tree_sequences_A.csv", optional: true), emit: rsv_a
    tuple path("tree_sequences_B.fasta", optional: true),
        path("tree_sequences_B.csv", optional: true), emit: rsv_b
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Generate the report fasta and csv files for generating the final report
    generate_csv_fasta.py \\
        --meta $manifest \\
        --version $version \\
        --output . \\
        --reference $reference_dir \\
        --coverage-threshold $igv_cutoff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
