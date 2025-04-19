//
// RSV GGene Genotyping Subworkflow
// This subworkflow performs genotyping of RSV samples using a series processes
//

include { BLAST_MAKEBLASTDB as GENOTYPING_MAKEBLASTDB } from '../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN as BLASTN_GENOTYPING           } from '../../modules/nf-core/blast/blastn/main'
include { CAT_CAT as CAT_FASTA                        } from '../../modules/nf-core/cat/cat/main'
include { MAFFT_ALIGN                                 } from '../../modules/nf-core/mafft/align/main'
include { FASTTREE                                    } from '../../modules/nf-core/fasttree/main'
include { EXTRACT_GENE_SEQUENCES                      } from '../../modules/local/extract_gene_sequences'
include { PARSE_GGENE_BLASTN                          } from '../../modules/local/parse_ggene_blastn'
include { VISUALIZE_PHYLOGENETIC_TREE                 } from '../../modules/local/visualize_tree/main'


workflow RSV_GGENE_GENOTYPING {

    take:
    ch_consensus_fasta         // Channel: [ val(meta), path(fasta) ]
    ch_genotyping_gg_blast_db  // Channel: [ val(meta), path(db)    ]
    ch_matched_ref_fasta       // Channel: [ val(meta), path(fasta) ]
    ch_matched_ref_gff         // Channel: [ val(meta), path(gff)   ]

    main:

    ch_versions = Channel.empty()

    //
    // Step 1: Extract g-genes from the assembled fasta file
    //
    EXTRACT_GENE_SEQUENCES (
        ch_consensus_fasta
            .join(ch_matched_ref_fasta, by: [0])
            .join(ch_matched_ref_gff, by: [0]),
        "CDS_7"
    )
    ch_ggene_consensus_fasta = EXTRACT_GENE_SEQUENCES.out.fasta
    ch_versions = ch_versions.mix(EXTRACT_GENE_SEQUENCES.out.versions.first())

    //
    // Step 2: Run BLAST against reference database for genotyping
    //
    // GENOTYPING_MAKEBLASTDB (
    //     ch_genotype_ggene_ref_fasta.map {
    //         meta, fasta ->
    //             [ meta + [id:"RSV_GGENE"], fasta ]
    //     }
    // )
    // ch_versions = ch_versions.mix(GENOTYPING_MAKEBLASTDB.out.versions)

    BLASTN_GENOTYPING (
        ch_ggene_consensus_fasta,
        ch_genotyping_gg_blast_db
        // GENOTYPING_MAKEBLASTDB.out.db
    )
    ch_blast_out = BLASTN_GENOTYPING.out.txt
    ch_versions  = ch_versions.mix(BLASTN_GENOTYPING.out.versions.first())

    //
    // Step 3: G-gene Genotyping
    //
    PARSE_GGENE_BLASTN ( ch_blast_out )
    ch_versions = ch_versions.mix(PARSE_GGENE_BLASTN.out.versions.first())

    emit:

    blast_out = ch_blast_out // tuple: [ val(meta), path(txt) ]
    ggene_consensus_fasta = ch_ggene_consensus_fasta  // tuple: [ val(meta), path(extracted_fasta) ]
    versions = ch_versions
}
