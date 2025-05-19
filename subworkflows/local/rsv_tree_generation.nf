//
// SUBWORKFLOW: RSV Tree Generation
//

include { CAT_CAT as CAT_FASTA        } from '../../modules/nf-core/cat/cat/main'
include { PHYLOGENY_TREE_GENERATION   } from '../../subworkflows/local/phylogeny_tree'


workflow RSV_TREE_GENERATION {

    take:
    ch_query_fasta    // Channel: [ val(meta), path(fasta) ]
    ch_tree_ref_files // Channel: [ val(subtype), val(out_grp), path(ref_fasta), path(ref_meta), path(color) ]

    main:

    ch_versions = Channel.empty()

    ch_query_fasta
        .map { meta, fasta ->
            [meta.subtype, meta, fasta]
        }
        .combine ( ch_tree_ref_files, by: [0] )
        .map { subtype, meta, query_fasta, out_grp, ref_fasta, ref_meta, color ->
            [ meta, query_fasta, out_grp, ref_fasta, ref_meta, color ]
        }
        .set { ch_phy_tree_data }

    ch_phy_tree_data
        .map {
            [ it[0], [it[1], it[3]] ] // Channel: [ meta, [query_fasta, ref_fasta] ]
    }.set{ ch_tree_fasta }

    //
    // MODULE: Concat query and reference fasta together
    //
    CAT_FASTA ( ch_tree_fasta )
    ch_concat_fasta = CAT_FASTA.out.file_out
    ch_versions     = ch_versions.mix(CAT_FASTA.out.versions.first())

    //
    // SUBWORKFLOW: Generate the phylogeny tree
    //
    ch_concat_fasta.join(
        ch_phy_tree_data.map {
            meta, query_fasta, out_grp, ref_fasta, ref_meta, color ->
                [ meta, out_grp, ref_meta, color ]
        },
        by: [0]
    )
    .set { ch_tree_input }

    PHYLOGENY_TREE_GENERATION ( ch_tree_input )
    ch_phy_tree_plot = PHYLOGENY_TREE_GENERATION.out.phy_tree_plot
    ch_versions = ch_versions.mix(PHYLOGENY_TREE_GENERATION.out.versions)

    emit:

    phy_tree_plot = ch_phy_tree_plot   // channel: [ val(meta), path(png) ]
    versions = ch_versions
}
