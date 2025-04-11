//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA                  } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF                    } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_META                   } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GENOTYPE_META          } from '../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_KMA_INDEX                } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_GISAID_DB                } from '../../modules/nf-core/untar/main'
include { KMA_INDEX                               } from '../../modules/nf-core/kma/index/main'
include { BLAST_MAKEBLASTDB as GISAID_MAKEBLASTDB } from '../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_MAKEBLASTDB as WHG_MAKEBLASTDB    } from '../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_MAKEBLASTDB as GGENE_MAKEBLASTDB  } from '../../modules/nf-core/blast/makeblastdb/main'


workflow PREPARE_REFERENCE_FILES {

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress RSV genome fasta file if required
    //
    if (params.rsv_fasta_pool) {
        if (params.rsv_fasta_pool.endsWith('.gz')) {
            GUNZIP_FASTA (
                [ [:], file(params.rsv_fasta_pool, type: 'file', checkIfExists: true) ]
            )
            ch_rsv_fasta_pool = GUNZIP_FASTA.out.gunzip.map {
                meta, fasta -> [ [id:fasta.Name], fasta ]
            }
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
            ch_rsv_fasta_pool = Channel.value(file(params.rsv_fasta_pool, type: 'file', checkIfExists: true)).map{ [[id:it.Name], it] }
        }
    } else {
        GUNZIP_FASTA (
            [ [:], file("${projectDir}/vendor/RefPool/RSV.fasta.gz") ]
        )
        ch_rsv_fasta_pool = GUNZIP_FASTA.out.gunzip.map {
            meta, fasta -> [ [id:fasta.Name], fasta ]
        }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    }

    //
    // Uncompress RSV genome gff file if required
    //
    if (params.rsv_gff_pool) {
        if (params.rsv_gff_pool.endsWith('.gz')) {
            GUNZIP_GFF (
                [ [:], file(params.rsv_gff_pool, type: 'file', checkIfExists: true) ]
            )
            ch_rsv_gff_pool = GUNZIP_GFF.out.gunzip
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_rsv_gff_pool = Channel.value(file(params.rsv_gff_pool, type: 'file', checkIfExists: true)).map{ [[:], it] }
        }
    } else {
        GUNZIP_GFF (
            [ [:], file("${projectDir}/vendor/RefPool/RSV.gff.gz") ]
        )
        ch_rsv_gff_pool = GUNZIP_GFF.out.gunzip
        ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
    }

    //
    // Uncompress RSV genome meta file if required
    //
    if (params.rsv_meta_pool) {
        if (params.rsv_meta_pool.endsWith('.gz')) {
            GUNZIP_META (
                [ [:], file(params.rsv_meta_pool, type: 'file', checkIfExists: true) ]
            )
            ch_rsv_meta_pool = GUNZIP_META.out.gunzip
            ch_versions = ch_versions.mix(GUNZIP_META.out.versions)
        } else {
            ch_rsv_meta_pool = Channel.value(file(params.rsv_meta_pool, type: 'file', checkIfExists: true)).map{ [[:], it] }
        }
    } else {
        GUNZIP_META (
            [ [:], file("${projectDir}/vendor/RefPool/RSV.csv.gz") ]
        )
        ch_rsv_meta_pool = GUNZIP_META.out.gunzip
        ch_versions = ch_versions.mix(GUNZIP_META.out.versions)
    }

    //
    // Prepare reference files required by KMA for identifying the best reference genome
    //
    ch_kma_index = Channel.empty()
    if (params.kma_index) {
        if (params.kma_index.endsWith('.tar.gz')) {
            UNTAR_KMA_INDEX (
                [ [:], file(params.kma_index, type: 'file') ]
            )
            ch_kma_index = UNTAR_KMA_INDEX.out.untar
            ch_versions  = ch_versions.mix(UNTAR_KMA_INDEX.out.versions)
        } else {
            ch_kma_index = Channel.value( file(params.kma_index, type: 'dir') )
                .map { it -> [ [id:it.Name], it ] }
        }
    } else {
        KMA_INDEX ( ch_rsv_fasta_pool )
        ch_kma_index = KMA_INDEX.out.index
        ch_versions  = ch_versions.mix(KMA_INDEX.out.versions)
    }

    //
    // Prepare the BLAST DB files for GISAID genome
    //
    ch_gisaid_blast_db = Channel.empty()
    if (params.blast_gisaid_db) {
        if (params.blast_gisaid_db.endsWith('.tar.gz')) {
            UNTAR_GISAID_DB (
                [ [:], file(params.blast_gisaid_db, type: 'file') ]
            )
            ch_gisaid_blast_db = UNTAR_GISAID_DB.out.untar
            ch_versions = ch_versions.mix(UNTAR_GISAID_DB.out.versions)
        } else {
            ch_gisaid_blast_db = Channel.value( file(params.gisaid_blast_db, type: 'dir') ).map{ [[:], it] }
        }
    } else {
        GISAID_MAKEBLASTDB (
            [ [id:"GISAID"], file("${projectDir}/vendor/GISAIDDB/GISAID.fasta.gz", type: 'file', checkIfExists: true) ]
        )
        ch_gisaid_blast_db = GISAID_MAKEBLASTDB.out.db
        ch_versions = ch_versions.mix(GISAID_MAKEBLASTDB.out.versions)
    }

    //
    // Prepare the reference(db) files for wholegenome genotyping
    //
    ch_genotype_whg_blast_db = Channel.empty()
    ch_genotype_whg_ref_meta = Channel.empty()
    if (!params.skip_genotyping && !params.skip_wholegenome_genotyping) {
        if (params.genotype_whole_genome_fasta) {
            ch_genotype_whg_ref_fasta = Channel.value(
                file(params.genotype_whole_genome_fasta, type: 'file', checkIfExists: true)
            ).map { [[:], it] }
        } else {
            ch_genotype_whg_ref_fasta = Channel.value(
                file("${projectDir}/vendor/genotype/NextStrain.fasta.gz", type: 'file', checkIfExists: true)
            ).map { [[:], it] }
        }

        WHG_MAKEBLASTDB (
            ch_genotype_whg_ref_fasta
                .map { meta, fasta -> [ meta + [id:"RSV_AB"], fasta ] }
        )
        ch_genotype_whg_blast_db = WHG_MAKEBLASTDB.out.db
        ch_versions = ch_versions.mix(WHG_MAKEBLASTDB.out.versions)

        if (params.genotype_whole_genome_meta) {
            if (params.genotype_whole_genome_meta.endsWith('.gz')) {
                GUNZIP_GENOTYPE_META (
                    [ [:], file(params.genotype_whole_genome_meta, type: 'file', checkIfExists: true) ]
                )
                ch_genotype_whg_ref_meta = GUNZIP_GENOTYPE_META.out.gunzip
                ch_versions = ch_versions.mix(GUNZIP_GENOTYPE_META.out.versions)
            } else {
                ch_genotype_whg_ref_meta = Channel.value(
                    file(params.genotype_whole_genome_meta, type: 'file', checkIfExists: true)
                ). map { [[:], it] }
            }
        } else {
            GUNZIP_GENOTYPE_META (
                [ [:], file("${projectDir}/vendor/genotype/NextStrain.tsv.gz", type: 'file', checkIfExists: true) ]
            )
            ch_genotype_whg_ref_meta = GUNZIP_GENOTYPE_META.out.gunzip
            ch_versions = ch_versions.mix(GUNZIP_GENOTYPE_META.out.versions)
        }
    }

    //
    // Prepare the reference(db) files for g-gene genotyping
    //
    ch_genotype_gg_blast_db = Channel.empty()
    if (!params.skip_genotyping && !params.skip_ggene_genotyping) {
        if (params.genotype_ggene_fasta) {
            ch_genotype_gg_ref_fasta = Channel.value(
                file(params.genotype_ggene_fasta, type: 'file', checkIfExists: true)
            ).map { [[:], it] }
        } else {
            ch_genotype_gg_ref_fasta = Channel.value(
                file("${projectDir}/vendor/genotype/G_subtype.fasta.gz", type: 'file', checkIfExists: true)
            ).map { [[:], it] }
        }

        GGENE_MAKEBLASTDB (
            ch_genotype_gg_ref_fasta.map {
                meta, fasta -> [ meta + [id:"RSV_GGENE"], fasta] }
        )
        ch_genotype_gg_blast_db = GGENE_MAKEBLASTDB.out.db
        ch_versions = ch_versions.mix(GGENE_MAKEBLASTDB.out.versions)
    }

    //
    // Prepare the whole genome tree reference files for generating phylogenetic tree
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
    .set { ch_tree_whg_ref }

    //
    // Prepare the g-gene tree reference files for generating phylogenetic tree
    //
    Channel.of(
        [
            "SubtypeA",
            "AY911262_USA_1956|GA1",
            file("${projectDir}/vendor/TreeReference/g_gene_representative_ref_A.fasta", type: 'file', checkIfExists: true),
            file("${projectDir}/vendor/TreeReference/g_gene_representative_ref_A.csv",   type: 'file', checkIfExists: true),
            file("${projectDir}/vendor/TreeReference/g_gene_color_A.csv",                type: 'file', checkIfExists: true)
        ],
        [
            "SubtypeB",
            "AY353550_USA_1977|GB1",
            file("${projectDir}/vendor/TreeReference/g_gene_representative_ref_B.fasta", type: 'file', checkIfExists: true),
            file("${projectDir}/vendor/TreeReference/g_gene_representative_ref_B.csv",   type: 'file', checkIfExists: true),
            file("${projectDir}/vendor/TreeReference/g_gene_color_B.csv",                type: 'file', checkIfExists: true)
        ]
    )
    .set { ch_tree_gg_ref }

    emit:
    rsv_fasta_pool          = ch_rsv_fasta_pool  // tuple: [ val(meta), path(fasta) ]
    rsv_gff_pool            = ch_rsv_gff_pool    // tuple: [ val(meta), path(gff)   ]
    rsv_meta_pool           = ch_rsv_meta_pool   // tuple: [ val(meta), path(csv)   ]
    kma_index               = ch_kma_index       // tuple: [ val(meta), path(index) ]
    gisaid_blast_db         = ch_gisaid_blast_db // tuple: [ val(meta), path(db)    ]

    // genotyping reference files
    genotype_whg_blast_db   = ch_genotype_whg_blast_db // tuple: [ val(meta), path(db)  ]
    genotype_whg_ref_meta   = ch_genotype_whg_ref_meta // tuple: [ val(meta), path(tsv) ]
    genotype_gg_blast_db    = ch_genotype_gg_blast_db  // tuple: [ val(meta), path(db)  ]

    // Phylogenetic tree generation files
    tree_whg_ref            = ch_tree_whg_ref // Channel: [ val(subtype), val(out_grp), path(ref_fasta), path(ref_meta), path(color) ]
    tree_gg_ref             = ch_tree_gg_ref  // Channel: [ val(subtype), val(out_grp), path(ref_fasta), path(ref_meta), path(color) ]

    versions                = ch_versions     // channel: [ versions.yml ]
}
