process VISUALIZE_PHYLOGENETIC_TREE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/stjudecab/rsvrecon:v4.3"

    input:
    tuple val(meta), path(tree), path(annotation), path(color)

    output:
    tuple val(meta), path("*.png"), emit: tree_plot
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in stjudecab/rsvrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp $annotation genotype.csv
    echo "${meta.id},Query,Query\\n" >> genotype.csv

    visualize_tree.R \\
        $tree \\
        "${meta.strain}" \\
        genotype.csv \\
        ${prefix}.phylogenetic.png \\
        $color \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-ggtree: \$(Rscript -e "library(ggtree); cat(as.character(packageVersion('ggtree')))")
        r-treeio: \$(Rscript -e "library(treeio); cat(as.character(packageVersion('treeio')))")
    END_VERSIONS
    """
}
