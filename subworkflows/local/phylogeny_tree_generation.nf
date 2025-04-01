//
// SUBWORKFLOW: Phylogeny Tree Generation
//

include { CAT_CAT as CAT_FASTA        } from '../../modules/nf-core/cat/cat/main'
include { MAFFT_ALIGN                 } from '../../modules/nf-core/mafft/align/main'
include { FASTTREE                    } from '../../modules/nf-core/fasttree/main'
include { VISUALIZE_PHYLOGENETIC_TREE } from '../../modules/local/visualize_tree/main'


workflow PHYLOGENY_TREE_GENERATION {

    take:
    ch_query_fasta    // Channel: [ val(meta), path(fasta) ]
    ch_tree_ref_files // Channel: [ val(subtype), val(out_grp), path(ref_fasta), path(ref_meta), path(color) ]

    main:

    ch_versions = Channel.empty()

    ch_query_fasta
        .map { meta, fasta ->
            [meta.subtype, meta, fasta]
        }
        .join ( ch_tree_ref_files, by: [0])
        .map { subtype, meta, query_fasta, out_grp, ref_fasta, ref_meta, color ->
            [ meta, fasta, out_grp, ref_fasta, ref_meta, color ]
        }
        .set { ch_phy_tree_input }

    ch_phy_tree_input
        .map {
            [ it[0], [it[1], it[3]] ] // Channel: [ meta, [query_fasta, ref_fasta] ]
    }.set{ ch_tree_fasta }

    //
    // MODULE: Concat query and reference fasta together
    //
    CAT_FASTA (ch_concat_fasta)
    ch_concat_fasta = CAT_FASTA.out.file_out
    ch_versions     = ch_versions.mix(CAT_FASTA.out.versions.first())

    //
    // MODULE: Run multiple nucleotide sequence alignment with MAFFT
    //
    MAFFT_ALIGN (
        ch_concat_fasta,
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
        .join( ch_phy_tree_input.map {
            meta, query_fasta, out_grp, ref_fasta, ref_meta, color ->
                [ meta, out_grp, ref_meta, color ]
            },
            by: [0]
        )
        .map {
            meta, tree, out_grp, ref_meta, color ->
                [ meta + [strain:out_grp], tree, ref_meta, color ]
        }
        .set { ch_tree_render_input }

    VISUALIZE_PHYLOGENETIC_TREE ( ch_tree_render_input )
    ch_versions = ch_versions.mix(VISUALIZE_PHYLOGENETIC_TREE.out.versions.first())

    emit:

    versions = ch_versions
}
