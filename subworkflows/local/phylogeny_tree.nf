//
// SUBWORKFLOW: Phylogeny Tree Generation
//

include { MAFFT_ALIGN                 } from '../../modules/nf-core/mafft/align/main'
include { FASTTREE                    } from '../../modules/nf-core/fasttree/main'
include { VISUALIZE_PHYLOGENETIC_TREE } from '../../modules/local/visualize_tree/main'


workflow PHYLOGENY_TREE_GENERATION {

    take:
    ch_tree_input  // channel: [ val(meta), path(fasta), val(outgrp), path(ref_meta), path(color) ]

    main:

    ch_versions = Channel.empty()

    ch_fasta = ch_tree_input.map {
        meta, fasta, outgrp, ref_meta, color ->
        [ meta, fasta ]
    }

    //
    // MODULE: Run multiple nucleotide sequence alignment with MAFFT
    //
    MAFFT_ALIGN (
        ch_fasta,
        [[:],[]],
        [[:],[]],
        [[:],[]],
        [[:],[]],
        [[:],[]],
        false
    )

    ch_mafft_fasta = MAFFT_ALIGN.out.fas
    ch_versions    = ch_versions.mix(MAFFT_ALIGN.out.versions.first())

    //
    // MODULE: Generate the Phylogeny tree in `NEWICK` Tree format with FastTree
    //
    FASTTREE ( ch_mafft_fasta )
    ch_phylogeny_tree = FASTTREE.out.phylogeny
    ch_versions = ch_versions.mix(FASTTREE.out.versions.first())

    //
    // Draw the phylogeny tree in PNG format
    //
    ch_phylogeny_tree
        .join(ch_tree_input.map {
            meta, fasta, outgrp, ref_meta, color ->
                [ meta, outgrp, ref_meta, color ]
            },
            by: [0]
        )
        .map {
            meta, tree, out_grp, ref_meta, color ->
                [ meta + [strain:out_grp], tree, ref_meta, color ]
        }
        .set { ch_tree_render_input }

    VISUALIZE_PHYLOGENETIC_TREE ( ch_tree_render_input )
    ch_phy_tree_plot = VISUALIZE_PHYLOGENETIC_TREE.out.tree_plot
    ch_versions = ch_versions.mix(VISUALIZE_PHYLOGENETIC_TREE.out.versions.first())

    emit:

    phy_tree_plot = ch_phy_tree_plot   // channel: [ val(meta), path(png) ]
    versions = ch_versions
}
