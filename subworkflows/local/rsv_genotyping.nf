//
// RSV Genotyping Subworkflow
// This subworkflow performs genotyping of RSV samples using a series processes
//

include { BLAST_MAKEBLASTDB as GENOTYPING_MAKEBLASTDB } from '../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN as BLASTN_GENOTYPING           } from '../../modules/nf-core/blast/blastn/main'
include { CAT_CAT as CAT_FASTA                        } from '../../modules/nf-core/cat/cat/main'
include { MAFFT_ALIGN                                 } from '../../modules/nf-core/mafft/align/main'
include { FASTTREE                                    } from '../../modules/nf-core/fasttree/main'
include { PARSE_BLASTN                                } from '../../modules/local/parse_blastn'
include { VISUALIZE_PHYLOGENETIC_TREE                 } from '../../modules/local/visualize_tree/main'


workflow RSV_GENOTYPING {

    take:
    ch_consensus_fasta       // Channel: [ val(meta), path(fasta) ]
    ch_genotyping_blast_ref  // path/to/RSV_AB.fasta
    ch_genotyping_blast_meta // path/to/RSV_AB_meta.txt

    main:

    ch_versions = Channel.empty()

    //
    // Step 1: Run BLAST against reference database for genotyping
    //
    GENOTYPING_MAKEBLASTDB (
        ch_genotyping_blast_ref.map {
            meta, fasta ->
                [ meta + [id:"RSV_AB"], fasta ]
        }
    )
    ch_versions = ch_versions.mix(GENOTYPING_MAKEBLASTDB.out.versions)

    BLASTN_GENOTYPING (
        ch_consensus_fasta,
        GENOTYPING_MAKEBLASTDB.out.db
    )
    ch_blast_out = BLASTN_GENOTYPING.out.txt
    ch_versions  = ch_versions.mix(BLASTN_GENOTYPING.out.versions.first())

    //
    // Step 2: Select appropriate
    //
    PARSE_BLASTN ( ch_blast_out, ch_genotyping_blast_meta.map{ it[1] } )
    ch_versions = ch_versions.mix(PARSE_BLASTN.out.versions.first())

    //
    // Step 3: Generate tree
    //

    Channel.of(
        [
            "SubtypeA",
            "MG642074|A",
            file("${projectDir}/vendor/TreeReference/representative_ref_A.fasta", type: 'file', checkIfExists: true),
            file("${projectDir}/vendor/TreeReference/representative_ref_A.csv",   type: 'file', checkIfExists: true),
            file("${projectDir}/vendor/TreeReference/color_A.csv",                type: 'file', checkIfExists: true)
        ],
        [
            "SubtypeB",
            "Ger/302/98-99|B",
            file("${projectDir}/vendor/TreeReference/representative_ref_B.fasta", type: 'file', checkIfExists: true),
            file("${projectDir}/vendor/TreeReference/representative_ref_B.csv",   type: 'file', checkIfExists: true),
            file("${projectDir}/vendor/TreeReference/color_B.csv",                type: 'file', checkIfExists: true)
        ]
    )
    .set { ch_tree_ref_files }

    ch_consensus_fasta
        .map { meta, fasta ->
            [meta.subtype, meta, fasta]
        }
        .join ( ch_tree_ref_files, by: [0])
        .map { subtype, meta, fasta, out_group, ref_fasta, ref_csv, color_csv ->
            [ meta, fasta, out_group, ref_fasta, ref_csv, color_csv ]
        }
        .set { ch_tree_input }

    ch_tree_input
        .map {
            [ it[0], [it[1], it[3]] ] // Channel: [ meta, [fasta, ref_fasta] ]
        }.set{ ch_concat_fasta }

    CAT_FASTA (ch_concat_fasta)
    ch_versions = ch_versions.mix(CAT_FASTA.out.versions.first())

    // Mafft to perform multiple sequence alignment
    MAFFT_ALIGN (
        CAT_FASTA.out.file_out,
        [[:],[]],
        [[:],[]],
        [[:],[]],
        [[:],[]],
        [[:],[]],
        false
    )

    ch_mafft_fas = MAFFT_ALIGN.out.fas
    ch_versions = ch_versions.mix(MAFFT_ALIGN.out.versions.first())

    // Fasttree to generate the tree file
    FASTTREE ( ch_mafft_fas )
    ch_phylogeny_tree = FASTTREE.out.phylogeny
    ch_versions = ch_versions.mix(FASTTREE.out.versions.first())


    // Join the tree reference files with the tree file
    ch_phylogeny_tree
        .join( ch_tree_input.map {
            meta, fasta, out_group, ref_fasta, ref_csv, color_csv ->
                [ meta, out_group, ref_csv, color_csv ]
            },
            by: [0]
        )
        .map {
            meta, tree, out_group, ref_csv, color_csv ->
                [ meta + [strain:out_group], tree, ref_csv, color_csv ]
        }
        .set { ch_render_tree }

    // Run the drawing script for generating trees
    VISUALIZE_PHYLOGENETIC_TREE ( ch_render_tree )
    ch_versions = ch_versions.mix(VISUALIZE_PHYLOGENETIC_TREE.out.versions.first())

    emit:

    versions = ch_versions
}
