#!/usr/bin/env Rscript

# Phylogenetic Tree Visualization Script
# This script generates a publication-ready phylogenetic tree visualization
# with strain annotations, clade colorization, and type indicators.

#Author: Haidong Yi (hyi@stjude.org)
#        Lei Li (lei.li@stjude.org)
#Date: March 28, 2025

# Load required libraries
suppressPackageStartupMessages({
    library(ggtree)
    library(treeio)
    library(ggplot2)
})

# read parameters
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
    stop("Usage: Rscript tree_visualization.R <tree_file> <outgroup_strain> <annotation_file> <output_figure> <color_file>")
}

tree_file <- args[1]
outgroup <- args[2]
anno_file <- args[3]
fig_name <- args[4]
color_file <- args[5]

# Log parameters for debugging
message("Parameters:")
message(paste("Tree file:", tree_file))
message(paste("Outgroup:", outgroup))
message(paste("Annotation file:", anno_file))
message(paste("Output figure:", fig_name))
message(paste("Color file:", color_file))

# Read input data
tryCatch({
    tree <- read.newick(tree_file)
    strain_data <- read.csv(anno_file, header = FALSE)
    color_scheme <- read.csv(color_file, header = FALSE)
}, error = function(e) {
    stop(paste("Error reading input files:", e$message))
})

# Process annotation data
num_strains <- nrow(strain_data)
colnames(strain_data) <- c('label', 'Clade', 'Type')
strain_data$Type <- factor(strain_data$Type, levels = c('Reference', 'Query'))

# Process color scheme
color_set <- color_scheme$V2
names(color_set) <- color_scheme$V1

# Reroot tree with specified outgroup
root_tree <- root(tree, outgroup = outgroup)

# Generate tree visualization
p <- ggtree(root_tree, size = 2) %<+% strain_data
node_data <- p$data
max_x <- max(node_data$x)

# Configure plot aesthetics
p <- p +
    geom_tiplab(size = 6, color = "black", geom = 'label') +
    geom_tippoint(aes(color = Clade, shape = Type, size = Type)) +
    geom_treescale(linesize = 1.5, color = 'black', x = max_x * 0.8) +
    theme(
        legend.position = 'right',
        legend.background = element_rect(),
        legend.key = element_blank(),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size = 15),
        title = element_text(size = 15)
    ) +
    xlim(0, max_x * 1.8) +
    scale_color_manual(values = color_set) +
    scale_shape_manual(values = c(20, 17)) +
    scale_size_manual(values = c(7, 9))

# Save plot with dynamic sizing based on strain count
message(paste("Saving figure with", num_strains, "strains"))

if (num_strains > 30) {
    size_factor <- sqrt(num_strains / 30)
    width <- 15 * size_factor
    height <- 15 * size_factor
    message(paste("Using scaled dimensions:", width, "x", height))
} else {
    width <- 15 * 1.2
    height <- 15 * 1.2
    message(paste("Using default dimensions:", width, "x", height))
}

ggsave(plot = p, filename = fig_name, width = width, height = height)
message(paste("Figure saved to", fig_name))

