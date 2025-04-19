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


workflow RSV_WHOLEGENOME_GENOTYPING {

    take:
    ch_consensus_fasta       // Channel: [ val(meta), path(fasta) ]
    ch_genotyping_whg_blast_db   // Channel: [ val(meta), path(db)    ]
    ch_genotyping_whg_blast_meta // Channel: [ val(meta), path(csv)   ]

    main:

    ch_versions = Channel.empty()

    //
    // Step 1: Run BLAST against reference database for genotyping
    //
    // GENOTYPING_MAKEBLASTDB (
    //     ch_genotyping_blast_ref.map {
    //         meta, fasta ->
    //             [ meta + [id:"RSV_AB"], fasta ]
    //     }
    // )
    // ch_versions = ch_versions.mix(GENOTYPING_MAKEBLASTDB.out.versions)

    BLASTN_GENOTYPING (
        ch_consensus_fasta,
        ch_genotyping_whg_blast_db
        // GENOTYPING_MAKEBLASTDB.out.db
    )
    ch_blast_out = BLASTN_GENOTYPING.out.txt
    ch_versions  = ch_versions.mix(BLASTN_GENOTYPING.out.versions.first())

    //
    // Step 2: Select appropriate
    //
    PARSE_BLASTN ( ch_blast_out, ch_genotyping_whg_blast_meta.map{ it[1] } )
    ch_versions = ch_versions.mix(PARSE_BLASTN.out.versions.first())

    emit:

    blast_out = ch_blast_out  // tuple: [ val(meta), path(txt) ]
    versions  = ch_versions
}
