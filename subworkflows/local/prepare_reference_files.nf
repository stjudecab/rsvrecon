//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA        } from '../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_KMA_INDEX      } from '../../modules/nf-core/untar/main'
include { CUSTOM_GETCHROMSIZES          } from '../../modules/nf-core/custom/getchromsizes/main'
include { KMA_INDEX                     } from '../../modules/nf-core/kma/index/main'



workflow PREPARE_REFERENCE_FILES {

    take:
    fasta  // from `param.fasta`

    main:

    ch_versions = Channel.empty()

    //
    // Umcompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA (
            [ [:], fasta ]
        )
        ch_fasta = GUNZIP_FASTA.out.gunzip.map {
            meta, fasta ->
                [ [id:fasta.Name], fasta ]
        }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(params.fasta)).map {
            meta, fasta ->
                [ [id:fasta.Name], fasta ]
        }
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai   // [ val(meta), path(fai) ]
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes // [ val(meta), path(sizes) ]
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Prepare reference files required by KMA for identifying best reference
    //
    ch_kma_index = Channel.empty()
    if (params.kma_index) {
        if (params.kma_index.endsWith('.tar.gz')) {
            UNTAR_KMA_INDEX (
                [ [:], params.kma_db ]
            )
            ch_kma_index = UNTAR_KMA_INDEX.out.untar
            ch_versions  = ch_versions.mix(UNTAR_KMA_INDEX.out.versions)
        } else {
            ch_kma_index = Channel.value(file(params.kma_index, type: 'dir', checkIfExists: true))
        }
    } else {
        KMA_INDEX ( ch_fasta )
        ch_kma_index = KMA_INDEX.out.index
        ch_versions  = ch_versions.mix(KMA_INDEX.out.versions)
    }

    emit:
    fasta             = ch_fasta       // tuple: [ val(meta), path(fasta) ]
    fai               = ch_fai         // tuple: [ val(meta), path(fai) ]
    chrom_sizes       = ch_chrom_sizes // tuple: [ val(meta), path(sizes) ]
    kma_index         = ch_kma_index   // tuple: [ val(meta), path("kmaindex/name") ]

    versions          = ch_versions    // channel: [ versions.yml ]
}
