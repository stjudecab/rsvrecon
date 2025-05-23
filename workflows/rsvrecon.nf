/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER CHECK AND PRESET
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check input path parameters
def checkPathParamList = [
    params.rsv_fasta_pool, params.rsv_gff_pool, params.rsv_meta_pool,
    params.kma_index, params.blast_gisaid_db,
    params.genotype_whole_genome_fasta, params.genotype_whole_genome_meta, params.genotype_ggene_fasta
]

for (param in checkPathParamList) {
    if (param) { file(param, checkIfExists: true) }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)   : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { KMA_MAP            } from '../modules/local/kma_map'
include { READ_KMA           } from '../modules/local/read_kma'
include { IGVTOOLS_COUNT     } from '../modules/local/igvtools_count'
include { ASSEMBLE_SEQUENCE  } from '../modules/local/assemble_sequence'
include { GENERATE_CSV_FASTA } from '../modules/local/generate_csv_fasta'
include { GENERATE_REPORT    } from '../modules/local/generate_report'

//
// SUBWORKFLOW: Consisting a mix of local and nf-core/modules
//
include { FASTQ_TRIM_FASTP_FASTQC    } from '../subworkflows/local/fastq_trim_fastp_fastqc'
include { PREPARE_REFERENCE_FILES    } from '../subworkflows/local/prepare_reference_files'
include { BAM_SORT_STATS_SAMTOOLS    } from '../subworkflows/local/bam_sort_stats_samtools/main'
include { RSV_WHOLEGENOME_GENOTYPING } from '../subworkflows/local/rsv_genotyping'
include { RSV_GGENE_GENOTYPING       } from '../subworkflows/local/rsv_ggene_genotyping'
include {
    RSV_TREE_GENERATION as PHY_WHG_TREE;
    RSV_TREE_GENERATION as PHY_GGENE_TREE } from '../subworkflows/local/rsv_tree_generation'
include {
    PHYLOGENY_TREE_GENERATION as PHY_RSV_A;
    PHYLOGENY_TREE_GENERATION as PHY_RSV_B  } from '../subworkflows/local/phylogeny_tree'

//
// MODULE: Installed directly from nf-core/modules (possibly with some patches)
//
include { CAT_FASTQ                    } from '../modules/nf-core/cat/fastq/main'
include { STAR_GENOMEGENERATE          } from '../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN                   } from '../modules/nf-core/star/align/main'
include { SAMTOOLS_MPILEUP             } from '../modules/nf-core/samtools/mpileup/main'
include { BLAST_BLASTN as BLAST_GISAID } from '../modules/nf-core/blast/blastn/main'
include { NEXTCLADE_DATASETGET         } from '../modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN                } from '../modules/nf-core/nextclade/run/main'
include { MULTIQC as MULTIQC_RAWQC     } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_FINALQC   } from '../modules/nf-core/multiqc/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rsvrecon_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RSVRECON {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    // init version, report and multiqc channel
    ch_versions = Channel.empty()
    ch_report_files = Channel.empty()
    ch_multiqc_report = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference data files used by the pipeline
    //
    PREPARE_REFERENCE_FILES ()

    ch_rsv_fasta_pool  = PREPARE_REFERENCE_FILES.out.rsv_fasta_pool
    ch_rsv_gff_pool    = PREPARE_REFERENCE_FILES.out.rsv_gff_pool
    ch_rsv_meta_pool   = PREPARE_REFERENCE_FILES.out.rsv_meta_pool
    ch_kma_index       = PREPARE_REFERENCE_FILES.out.kma_index
    ch_gisaid_blast_db = PREPARE_REFERENCE_FILES.out.gisaid_blast_db

    ch_versions = ch_versions.mix( PREPARE_REFERENCE_FILES.out.versions )

    // group multiple fastq files from the same sample
    ch_samplesheet
        .map {
            meta, fastq ->
                meta.id = meta.id.split('_')[0..-2].join('_')
                [ meta, fastq ]
        }
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }  // ch_fastq of two branches: single (one replicate) and multiple (multiple replicate) )

    //
    // MODULE: Concatenate FastQ files from the same sample if required
    //
    CAT_FASTQ ( ch_fastq.multiple )
        .reads
        .mix ( ch_fastq.single )
        .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQ_TRIM_FASTP_FASTQC (
        ch_cat_fastq,
        [],
        false,
        params.save_trimmed_fail,
        false
    )
    ch_trimmed_json = FASTQ_TRIM_FASTP_FASTQC.out.trim_json
    ch_trimmed_fastq = FASTQ_TRIM_FASTP_FASTQC.out.reads
    ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)

    // multiqc files from fastqc and fastp
    ch_fastp_multiqc          = FASTQ_TRIM_FASTP_FASTQC.out.fastp_multiqc
    ch_fastqc_multiqc_pretrim = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_multiqc_pretrim
    ch_fastqc_multiqc_postrim = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_multiqc_postrim
    ch_fail_reads_multiqc     = FASTQ_TRIM_FASTP_FASTQC.out.fail_reads_multiqc

    //
    // MODULE: Run KMA for identifying best matched RSV genome from the pool
    //
    KMA_MAP ( ch_trimmed_fastq, ch_kma_index.map{ it[1] } )
    ch_kma_map_res = KMA_MAP.out.res
    ch_versions = ch_versions.mix(KMA_MAP.out.versions.first())

    //
    // READ the KMA results and get the best matched reference id
    //
    READ_KMA (
        ch_kma_map_res,
        ch_rsv_fasta_pool.map { it[1] },
        ch_rsv_gff_pool.map { it[1] },
        ch_rsv_meta_pool.map { it[1] }
    )
    ch_versions = ch_versions.mix(READ_KMA.out.versions.first())

    // move env vars to meta
    READ_KMA.out.fasta
        .map { meta, fasta, ref_subtype_env ->
            [meta + [subtype: ref_subtype_env], fasta]
        }
        .set { ch_matched_ref_fasta }

    READ_KMA.out.gff
        .map { meta, gff, ref_subtype_env ->
            [meta + [subtype: ref_subtype_env], gff]
        }
        .set { ch_matched_ref_gff }

    ch_trimmed_fastq.join(
        READ_KMA.out.fasta, by: [0]
    ).map { meta, fastq, fasta, ref_subtype_env ->
        [ meta + [subtype: ref_subtype_env], fastq ]
    }.set { ch_trimmed_fastq }

    //
    // SUBWORKFLOW: STAR index and mapping the matched genome
    //
    STAR_GENOMEGENERATE ( ch_matched_ref_fasta, [[:],[]] )
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())

    ch_star_multiqc = Channel.empty()
    STAR_ALIGN (
        ch_trimmed_fastq
            .join(STAR_GENOMEGENERATE.out.index, by: [0]),
        [[:], []],
        true,
        [],
        []
    )
    ch_star_bam = STAR_ALIGN.out.bam
    ch_star_multiqc = ch_star_multiqc.mix(STAR_ALIGN.out.log_final)
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    //
    // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
    //
    BAM_SORT_STATS_SAMTOOLS ( ch_star_bam, ch_matched_ref_fasta )

    ch_star_sorted_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
    ch_star_sorted_bai = BAM_SORT_STATS_SAMTOOLS.out.bai

    // these stats go for multiqc
    ch_star_sorted_stats    = BAM_SORT_STATS_SAMTOOLS.out.stats
    ch_star_sorted_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
    ch_star_sorted_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    //
    // MODULE: SAMTOOLS MPILEUP
    //
    SAMTOOLS_MPILEUP ( ch_star_sorted_bam.join(ch_matched_ref_fasta, by: [0]) )
    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions.first())

    //
    // MODULE: IGV-Tools to count
    //
    IGVTOOLS_COUNT ( ch_star_sorted_bam.join(ch_matched_ref_fasta, by: [0]) )
    ch_base_count = IGVTOOLS_COUNT.out.count
    ch_versions = ch_versions.mix(IGVTOOLS_COUNT.out.versions.first())

    //
    // MODULE: Assemble sequences
    //
    ASSEMBLE_SEQUENCE (
        ch_base_count.join(ch_matched_ref_fasta, by: [0]), params.igv_cutoff
    )
    ch_consensus_fasta = ASSEMBLE_SEQUENCE.out.consensus
    ch_versions = ch_versions.mix(ASSEMBLE_SEQUENCE.out.versions.first())

    //
    // MODULE: BLAST against the GISAID Database
    //
    BLAST_GISAID ( ch_consensus_fasta, ch_gisaid_blast_db )
    ch_versions = ch_versions.mix(BLAST_GISAID.out.versions.first())

    //
    // MODULE: NextClade to call the variants
    //

    // Get the NextClade RSV dataset based on ref_id
    NEXTCLADE_DATASETGET( ch_consensus_fasta.map{ it[0] }, [] )
    ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions.first())

    // Run the Nextclade to call variants
    NEXTCLADE_RUN (
        ch_consensus_fasta.join(
            NEXTCLADE_DATASETGET.out.dataset,
            by: [0]
        )
    )
    ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())

    //
    // SUBWORKFLOW: Genotyping the assembled genome
    //
    if (!params.skip_genotyping) {

        // MODULE: Genotyping with the whole genome
        if (!params.skip_wholegenome_genotyping) {

            ch_genotype_whg_blast_db = PREPARE_REFERENCE_FILES.out.genotype_whg_blast_db
            ch_genotype_whg_ref_meta = PREPARE_REFERENCE_FILES.out.genotype_whg_ref_meta

            // MODULE: Genotyping with the whole genome
            RSV_WHOLEGENOME_GENOTYPING (
                ch_consensus_fasta,
                ch_genotype_whg_blast_db,
                ch_genotype_whg_ref_meta
            )
            ch_versions = ch_versions.mix(RSV_WHOLEGENOME_GENOTYPING.out.versions)

            // Generate the whole genome phylogenetic tree
            ch_tree_whg_ref = PREPARE_REFERENCE_FILES.out.tree_whg_ref

            PHY_WHG_TREE ( ch_consensus_fasta, ch_tree_whg_ref )
            ch_versions = ch_versions.mix(PHY_WHG_TREE.out.versions)
        }

        if (!params.skip_ggene_genotyping) {

            ch_genotype_gg_blast_db = PREPARE_REFERENCE_FILES.out.genotype_gg_blast_db

            // MODULE: Genotyping with the G-gene only
            RSV_GGENE_GENOTYPING (
                ch_consensus_fasta,
                ch_genotype_gg_blast_db,
                ch_matched_ref_fasta,
                ch_matched_ref_gff
            )
            ch_versions = ch_versions.mix(RSV_GGENE_GENOTYPING.out.versions)

            // Generate the ggene phylogenetic tree
            ch_tree_gg_ref = PREPARE_REFERENCE_FILES.out.tree_gg_ref

            PHY_GGENE_TREE (
                RSV_GGENE_GENOTYPING.out.ggene_consensus_fasta, ch_tree_gg_ref
            )
            ch_versions = ch_versions.mix(PHY_GGENE_TREE.out.versions)
        }
    }

    if (!params.skip_report && !params.skip_qc && !params.skip_fastp && !params.skip_genotyping
        && !params.skip_wholegenome_genotyping && !params.skip_ggene_genotyping)
    {

        ch_report_files = ch_star_sorted_flagstat
            .join(ch_consensus_fasta                           , by: [0])
            .join(ch_matched_ref_fasta                         , by: [0])
            .join(ch_matched_ref_gff                           , by: [0])
            .join(ch_base_count                                , by: [0])
            .join(NEXTCLADE_RUN.out.csv                        , by: [0])
            .join(RSV_WHOLEGENOME_GENOTYPING.out.whg_genotype  , by: [0])
            .join(RSV_WHOLEGENOME_GENOTYPING.out.blast_out     , by: [0])
            .join(RSV_GGENE_GENOTYPING.out.gg_genotype         , by: [0])
            .join(RSV_GGENE_GENOTYPING.out.blast_out           , by: [0])
            .join(BLAST_GISAID.out.txt                         , by: [0])
            .join(
                PHY_WHG_TREE
                    .out
                    .phy_tree_plot
                    .map { meta, png ->
                        def new_meta = meta.clone()
                        new_meta.remove('strain')
                        [new_meta, png]
                    },
                by: [0]
            )
            .join(
                PHY_GGENE_TREE
                    .out
                    .phy_tree_plot
                    .map { meta, png ->
                        def new_meta = meta.clone()
                        new_meta.remove('strain')
                        [new_meta, png]
                    },
                by: [0]
            )

        // ch_report_files.view()

        ch_report_files = ch_trimmed_json
            .join(ch_kma_map_res, by: [0])
            .map { meta, json, res -> [meta.id, json, res] }
            .join(
                ch_report_files.map{
                    def sample_id = it[0].id
                    def files = it[1..-1]
                    return [sample_id] + files
                },
                by: [0]
            )

        // generate the manifest file
        ch_report_files.map {
            def sample_id = it[0]
            def files = it[1..-1]
            def manifest = file("${workDir}/tmp/${sample_id}_manifest.tsv")
            header = ['sample_id', 'qc_fastp', 'kma_out', 'flagstat', 'assembly', 'ref_fasta', 'ref_gff', 'igv_out', 'nextclade_out', 'whg_genotype', 'whg_blastout', 'ggene_genotype', 'ggene_blastout', 'blast_gisaid', 'whg_figure', 'ggene_figure']
            def file_paths = []
            files.each { file ->
                file_paths << "${file.toAbsolutePath()}"
            }
            manifest.text = "${header.join("\t")}\n"
            manifest.append("${sample_id}\t${file_paths.join('\t')}\n")
            return manifest
        }
        .collectFile (name: "manifest.tsv", sort: true, keepHeader: true)
        .set { ch_manifest_file }

        // get the workflow version
        Channel.of("${getWorkflowVersion()}")
            .collectFile(name: "version.txt", newLine: true)
            .set { ch_workflow_version }

        // Run the report generation module
        GENERATE_CSV_FASTA (
            ch_manifest_file,
            ch_workflow_version,
            file("${projectDir}/vendor", type: 'dir', checkIfExists: true),
            params.igv_cutoff
        )
        ch_versions = ch_versions.mix(GENERATE_CSV_FASTA.out.versions)

        // Generate the phylogenetic tree of RSV A&B
        PHY_RSV_A (
            GENERATE_CSV_FASTA.out.rsv_a.map {
                fasta, csv ->
                    [[id:RSV_A], fasta, "MG642074|A", csv, file("${projectDir}/vendor/TreeReference/color_A.csv")]
            }
        )
        ch_versions = ch_versions.mix(PHY_RSV_A.out.versions)

        PHY_RSV_B (
            GENERATE_CSV_FASTA.out.rsv_b.map {
                fasta, csv ->
                    [[id:RSV_B], fasta, "Ger/302/98-99|B", csv, file("${projectDir}/vendor/TreeReference/color_B.csv")]
            }
        )
        ch_versions = ch_versions.mix(PHY_RSV_B.out.versions)

        // Generate the final report
        GENERATE_REPORT (
            GENERATE_CSV_FASTA.out.report,
            file("${projectDir}/vendor", type: 'dir', checkIfExists: true),
            file("${projectDir}/assets", type: 'dir', checkIfExists: true),
            ch_manifest_file,
            file("${projectDir}/assets/report_logo.png", type: 'file', checkIfExists: true),
            ch_workflow_version,
            PHY_RSV_A.out.phy_tree_plot.map{it[1]}.ifEmpty([]),
            PHY_RSV_B.out.phy_tree_plot.map{it[1]}.ifEmpty([]),
            params.igv_cutoff
        )
        ch_versions = ch_versions.mix(GENERATE_REPORT.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    if (!params.skip_qc && !params.skip_multiqc) {

        //
        // MODULE: MultiQC for raw data
        //
        ch_methods_description      = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_methods_description_file = ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true)

        ch_multiqc_rawqc_files = Channel.empty()

        // pre-trim qc files
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_fastqc_multiqc_pretrim.collect().ifEmpty([]))

        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_collated_versions)
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_methods_description_file)

        MULTIQC_RAWQC (
            ch_multiqc_rawqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )

        //
        // MODULE: MultiQC for final pipeline outputs
        //
        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

        ch_multiqc_finalqc_files = Channel.empty()
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_collated_versions)
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_methods_description_file)

        // post-trim qc files
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postrim.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fail_reads_multiqc.collectFile(name: 'fail_mapped_reads_mqc.tsv').ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastp_multiqc.collect().ifEmpty([]))

        // post-alignment qc files
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_star_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_star_sorted_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_star_sorted_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_star_sorted_idxstats.collect{it[1]}.ifEmpty([]))

        MULTIQC_FINALQC (
            ch_multiqc_finalqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )

        ch_multiqc_report = MULTIQC_FINALQC.out.report
        ch_versions = ch_versions.mix(MULTIQC_FINALQC.out.versions)
    }

    emit:
    multiqc_report = ch_multiqc_report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
