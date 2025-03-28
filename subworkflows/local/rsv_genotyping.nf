//
// RSV Genotyping Subworkflow
// This subworkflow performs genotyping of RSV samples using a series processes
//

include { BLAST_MAKEBLASTDB as GENOTYPING_MAKEBLASTDB } from '../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN as BLASTN_GENOTYPING           } from '../../modules/nf-core/blast/blastn/main'
include { PARSE_BLASTN                                } from '../../modules/local/parse_blastn'


workflow RSV_GENOTYPING {

    take:
    ch_consensus_fasta       // Channel: [ val(meta), path(fasta) ]
    ch_genotyping_blast_ref  // path/to/RSV_AB.fasta
    ch_genotyping_blast_meta // path/to/RSV_AB_meta.txt

    main:

    ch_versions = Channel.empty()

    // Step 1: Run BLAST against reference database for genotyping
    BLAST_MAKEBLASTDB (
        ch_genotyping_blast_ref.map {
            meta, fasta ->
                [ meta + [id:"RSV_AB"], fasta ]
        }
    )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    BLAST_GENOTYPING (
        ch_consensus_fasta,
        BLAST_MAKEBLASTDB.out.db
    )
    ch_blast_out = BLAST_GENOTYPING.out.txt
    ch_versions = ch_versions.mix(BLAST_GENOTYPING.out.versions.first())

    // Step 2: Select appropriate
    PARSE_BLASTN (
        ch_blast_out.combine(
            ch_genotyping_blast_meta.map { it[1] }
        )
    )

    emit:

    versions = ch_versions
}
