process GENERATE_REPORT {
    tag "GENERATE REPORT"
    label 'process_low'

    conda "conda-forge::python=3.9.12 conda-forge::biopython=1.79 conda-forge::pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' :
        'quay.io/biocontainers/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' }"

    input:
    path("manifest/*")
    path(ref_candidate_dir)
    path(tree_ref_dir)
    path(f_mutation_dir)

    output:
    path "report.html"   , optional: true, emit: report_html
    path "report.pdf"    , optional: true, emit: report_pdf
    path "versions.yml"  , emit: versions

    script:
    """
    # Generate a unifed report using the staged folder and manifest
    # generate_report.py \\
    #     --workdir ./ \\
    #     --manifest manifest.tsv \\
    #     --output_html report.html \\
    #     --output_pdf report.pdf
    python --help

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
