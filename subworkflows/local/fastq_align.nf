//
// FastQ Alignment
//

include { BWA_INDEX           } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM             } from '../../modules/nf-core/bwa/mem/main'
include { STAR_GENOMEGENERATE } from '../../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN          } from '../../modules/nf-core/star/align/main'


workflow FASTQ_ALIGNMENT {

    take:
    ch_reads      // channel: [ val(meta), [ path(reads) ] ]
    ch_fasta      // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()
    ch_fastq_align_multiqc = Channel.empty()

    if (params.aligner == 'bwa') {
        //
        // Map reads with BWA
        //
        BWA_INDEX ( ch_fasta )
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

        BWA_MEM (
            ch_reads.join(BWA_INDEX.out.index, by: [0]),
            [[:], []],
            false
        )

        ch_bam = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    } else if (params.aligner == 'star') {
        //
        // Map reads with STAR
        //
        STAR_GENOMEGENERATE ( ch_fasta, [[:],[]] )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())

        STAR_ALIGN (
            ch_reads
                .join(STAR_GENOMEGENERATE.out.index, by: [0]),
            [[:], []],
            true,
            [],
            []
        )

        ch_bam = STAR_ALIGN.out.bam
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

        ch_fastq_align_multiqc = ch_fastq_align_multiqc.mix(STAR_ALIGN.out.log_final)

    } else {
        log.error "RSVrecon pipeline only supports aligner: `bwa` or `star`, but got ${params.aligner}."
    }

    emit:
    bam                 = ch_bam                 // channel: [ val(meta), path(bam) ]
    fastq_align_multiqc = ch_fastq_align_multiqc // channel: [ val(meta), path(log) ]

    versions            = ch_versions            // channel: [ path(versions.yml) ]
}
