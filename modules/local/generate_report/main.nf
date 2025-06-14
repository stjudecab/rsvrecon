process GENERATE_REPORT {
    tag "GENERATE_REPORT"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/stjudecab/rsvrecon_report:v0.1"

    input:
    path(report_csv)
    path(reference_dir)
    path(asset_dir)
    path(meta, stageAs: "manifest/*")
    path(logo, stageAs: "manifest/*")
    path(version, stageAs: "manifest/*")
    path(tree_a, stageAs: "trees/*")
    path(tree_b, stageAs: "trees/*")
    val igv_cutoff

    output:
    path "*.pdf"       , emit: report_pdf
    path "*.html"      , emit: report_html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def tree_a_cmd = tree_a ? "--tree-a $tree_a" : ''
    def tree_b_cmd = tree_b ? "--tree-b $tree_b" : ''
    """
    # create tmp folder
    mkdir -p temp

    # generate final report pdf and html file
    generate_report.py \\
        --csv $report_csv \\
        --meta $meta \\
        --output-dir . \\
        --temp-dir ./temp \\
        --asset-dir $asset_dir \\
        --version $version \\
        --logo $logo \\
        $tree_a_cmd \\
        $tree_b_cmd \\
        --igv-cutoff $igv_cutoff \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
