/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER CHECK AND PRESET
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*

// Check input path parameters
def checkPathParamList = [
    params.fasta, params.kma_index
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
include { KMA_MAP                    } from '../modules/local/kma_map'

//
// SUBWORKFLOW: Consisting a mix of local and nf-core/modules
//
include { FASTQ_TRIM_FASTP_FASTQC    } from '../subworkflows/local/fastq_trim_fastp_fastqc'


//
// MODULE: Installed directly from nf-core/modules (possibly with some patches)
//
include { CAT_FASTQ                  } from '../modules/nf-core/cat/fastq/main'
include { MULTIQC as MULTIQC_RAWQC   } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_FINALQC } from '../modules/nf-core/multiqc/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
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

    // init version and multiqc channel
    ch_versions = Channel.empty()
    ch_multiqc_report = Channel.empty()

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
    ch_trimmed_fastq = FASTQ_TRIM_FASTP_FASTQC.out.reads
    ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)

    // multiqc files from fastqc and fastp
    ch_fastp_multiqc          = FASTQ_TRIM_FASTP_FASTQC.out.fastp_multiqc
    ch_fastqc_multiqc_pretrim = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_multiqc_pretrim
    ch_fastqc_multiqc_postrim = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_multiqc_postrim
    ch_fail_reads_multiqc     = FASTQ_TRIM_FASTP_FASTQC.out.fail_reads_multiqc

    //
    // SUBWORKFLOW: Prepare reference genomes
    //
    PREPARE_REFERENCE_FILES ( param.fasta )

    fasta = PREPARE_REFERENCE_FILES.out.fasta
    kma_index = PREPARE_REFERENCE_FILES.out.kma_index

    ch_versions = ch_versions.mix( PREPARE_REFERENCE_FILES.out.versions )

    //
    // MODULE: Alignment with KMA for identifying best
    //
    KMA_MAP ( ch_trimmed_fastq, kma_index.map{ it[1] } )
    ch_versions = ch_versions.mix(KMA_MAP.out.versions.first())

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

        MULTIQC_RAWQC(
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
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_methods_description_file)

        // post-trim qc files
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postrim.collect().ifEmpty([]))

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
