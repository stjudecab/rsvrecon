//
// Read QC and trimming
//

include { FASTQC as FASTQC_RAW  } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../modules/nf-core/fastqc/main'
include { FASTP                 } from '../../modules/nf-core/fastp/main'

//
// Function that parses fastp json output to get total number of reads after trimming
//

import groovy.json.JsonSlurper

def getFastpReadsAfterFiltering(json_file) {
    def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toInteger()
}

def getFastpReadsBeforeFiltering(json_file) {
    def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
    return json['before_filtering']['total_reads'].toInteger()
}

def multiqcTsvFromList(tsv_data, header) {
    def tsv_string = ""
    if (tsv_data.size() > 0) {
        tsv_string += "${header.join('\t')}\n"
        tsv_string += tsv_data.join('\n')
    }
    return tsv_string
}


def fail_mapped_reads = [:]

workflow FASTQ_TRIM_FASTP_FASTQC {

    take:
    reads                // channel: [ val(meta), [reads] ]
    adapter_fasta        // file: adapter.fasta
    discard_trimmed_pass // value: boolean
    save_trimmed_fail    // value: boolean
    save_merged          // value: boolean

    main:

    ch_versions = Channel.empty()

    fastqc_raw_html = Channel.empty()
    fastqc_raw_zip  = Channel.empty()
    fastqc_multiqc_pretrim = Channel.empty()
    fastqc_multiqc_postrim = Channel.empty()
    if (!params.skip_qc && !params.skip_fastqc) {
        FASTQC_RAW ( reads )

        fastqc_raw_html = FASTQC_RAW.out.html
        fastqc_raw_zip  = FASTQC_RAW.out.zip

        fastqc_raw_zip
            .map { it -> it[1] }
            .set { fastqc_raw_zip_only }
        fastqc_raw_html
            .map { it -> it[1] }
            .set { fastqc_raw_html_only }
        fastqc_multiqc_pretrim = fastqc_multiqc_pretrim.mix( fastqc_raw_zip_only, fastqc_raw_html_only )
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first().ifEmpty(null))
    }

    trim_reads         = reads
    trim_json          = Channel.empty()
    trim_html          = Channel.empty()
    trim_log           = Channel.empty()
    trim_reads_fail    = Channel.empty()
    trim_reads_merged  = Channel.empty()
    fastqc_trim_html   = Channel.empty()
    fastqc_trim_zip    = Channel.empty()
    fastp_multiqc      = Channel.empty()
    fail_reads_multiqc = Channel.empty()

    if (!params.skip_fastp) {
        FASTP (
            reads,
            adapter_fasta,
            discard_trimmed_pass,
            save_trimmed_fail,
            save_merged)

        trim_reads        = FASTP.out.reads
        trim_json         = FASTP.out.json
        trim_html         = FASTP.out.html
        trim_log          = FASTP.out.log
        trim_reads_fail   = FASTP.out.reads_fail
        trim_reads_merged = FASTP.out.reads_merged

        fastp_multiqc     = trim_json.map { it -> it[1] }
        ch_versions       = ch_versions.mix(FASTP.out.versions.first().ifEmpty(null))

        //
        // Filter empty FastQ files after adapter trimming so FastQC doesn't fail
        //
        trim_reads
            .join(trim_json)
            .map {
                meta, reads, json ->
                    pass = getFastpReadsAfterFiltering(json) > 0
                    [ meta, reads, json, pass ]
            }
            .set { trim_reads_pass_fail }

        trim_reads_pass_fail
            .map { meta, reads, json, pass -> if (pass) [ meta, reads ] }
            .set { trim_reads }

        trim_reads_pass_fail
            .map {
                meta, reads, json, pass ->
                if (!pass) {
                    fail_mapped_reads[meta.id] = 0
                    num_reads = getFastpReadsBeforeFiltering(json)
                    return [ "$meta.id\t$num_reads" ]
                }
            }
            .collect()
            .map {
                tsv_data ->
                    def header = ['Sample', "Reads before trimming"]
                    multiqcTsvFromList(tsv_data, header)
            }
            .set { fail_reads_multiqc }

        if (!params.skip_qc && !params.skip_fastqc) {
            FASTQC_TRIM ( trim_reads )

            fastqc_trim_zip  = FASTQC_TRIM.out.zip
            fastqc_trim_html = FASTQC_TRIM.out.html
            fastqc_trim_zip
                .map { it -> it[1] }
                .set { fastqc_trim_zip_only }
            fastqc_trim_html
                .map { it -> it[1] }
                .set { fastqc_trim_html_only }
            fastqc_multiqc_postrim = fastqc_multiqc_postrim.mix( fastqc_trim_zip_only, fastqc_trim_html_only )
            ch_versions = ch_versions.mix(FASTQC_TRIM.out.versions.first().ifEmpty(null))
        }
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]
    trim_json          // channel: [ val(meta), [ json ] ]
    trim_html          // channel: [ val(meta), [ html ] ]
    trim_log           // channel: [ val(meta), [ log ] ]
    trim_reads_fail    // channel: [ val(meta), [ fastq.gz ] ]
    trim_reads_merged  // channel: [ val(meta), [ fastq.gz ] ]
    fail_reads_multiqc // val(string): tab-separated table string
    fastp_multiqc      // channel: [ [ json ] ]

    fastqc_raw_html    // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip     // channel: [ val(meta), [ zip ] ]
    fastqc_trim_html   // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip    // channel: [ val(meta), [ zip ] ]
    fastqc_multiqc_pretrim
    fastqc_multiqc_postrim

    versions = ch_versions // channel: [ versions.yml ]
}
