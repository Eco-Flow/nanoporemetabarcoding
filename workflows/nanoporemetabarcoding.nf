/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                         } from '../modules/nf-core/fastqc/main'
include { MULTIQC                        } from '../modules/nf-core/multiqc/main'
include { PORECHOP_PORECHOP              } from '../modules/nf-core/porechop/porechop/main'
include { NANOPLOT                       } from '../modules/nf-core/nanoplot/main'
include { NANOFILT                       } from '../modules/nf-core/nanofilt/main'
include { SEQKIT_SEQ as SEQKIT_REVCOMP_A } from '../modules/nf-core/seqkit/seq/main'
include { SEQKIT_SEQ as SEQKIT_REVCOMP_B } from '../modules/nf-core/seqkit/seq/main'
include { CUTADAPT as CUTADAPT_F         } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_R         } from '../modules/nf-core/cutadapt/main'
include { paramsSummaryMap               } from 'plugin/nf-schema'
include { paramsSummaryMultiqc           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText         } from '../subworkflows/local/utils_nfcore_nanoporemetabarcoding_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOPOREMETABARCODING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_input = ch_samplesheet
             | map {
                meta, fastq ->
                [[id:meta.id, single_end:true], fastq] // Needs to be declared as single end to run PYCHOPPER
             }

    //
    // MODULE: Run Nanoplot
    //

    // Check modules.config for arguments to pass to Nanofilt

    NANOFILT (
        ch_input,
        []
    )

    //
    // MODULE: Run CUTADAPT
    //

    // Run cutadapt to demultiplex and trim reads based on forward barcodes (tag + primer)
    // set in the nextflow.config file

    CUTADAPT_F (
        ch_input
    )

    // After trimmed reads are demultiplexed and barcodes are removed, reads are concatenated
    // into a single FASTQ file. FASTQ files are then reverse complemented using seqkit
    // and cutadapt is run again to trim and demultiplex based on reverse barcodes

    //ch_reads = CUTADAPT.out.reads
    //         | map {
    //            meta, fastqs ->
    //            fastqs
    //         }
    //         | flatten()
    //         | collectFile(name: "all_reads.fastq.gz", storeDir: "${workDir}/tmp")
    //         | map { file -> [[id:'all_reads', single_end:true], file] }


    // Flatten the output channel (FASTQs) from cutadapt demultiplex into indidual channels (FASTQ)
    // (check for the function flattenAndMap in the functions section)
    ch_input_f = flattenAndMap(CUTADAPT_F.out.reads)

    //
    // MODULE: Run SEQKIT reverse complement
    //

    // Barcodes are attached to both ends of the reads, so we need to reverse complement the reads
    // to trim and demultiplex based the other end

    SEQKIT_REVCOMP_A (
        ch_input_f
    )

    // Run cutadapt on the reverse complemented reads

    CUTADAPT_R (
        SEQKIT_REVCOMP_A.out.fastx
    )

    // Filter out FASTQs with less than 10 reads
    ch_input_filte = CUTADAPT_R.out.reads
                   | filter { meta, fastq ->
                        def count = fastq.countFastq()
                        count > 10 // Filter out FASTQs with less than 1000 reads
                   }


    // Reverse complement the reads again to get back to the original orientation
    SEQKIT_REVCOMP_B (
        ch_input_f
    )

    // Prepare raw, cleaned and demultiplexed reads for Nanoplot

    ch_raw       = ch_input
                 | map {
                    meta, fastq ->
                    [[id:"raw_${meta.id}"], fastq]
                 }

    ch_filt     = NANOFILT.out.filtreads
                 | map {
                    meta, fastq ->
                    [[id:"filt_${meta.id}"], fastq]
                 }

    //
    // MODULE: Run Nanoplot
    //

    NANOPLOT (
        ch_raw.mix(ch_filt).mix(CUTADAPT_F.out.reads)
    )

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_input
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'nanoporemetabarcoding_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    //MULTIQC (
    //    ch_multiqc_files.collect(),
    //    ch_multiqc_config.toList(),
    //    ch_multiqc_custom_config.toList(),
    //    ch_multiqc_logo.toList(),
    //    [],
    //    []
    //)

    emit:
    //multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// When demultiplexing, cutadapt emits a channel with FASTQs, but the next modules
// input single values. Use flattenAndMap so that each FASTQ is emitted seprately
// Function to flatten output channel (FASTQs) from cutadapt demultiplex into indidual channels (FASTQ)

def flattenAndMap(ch_fastqs) {
    ch_fastq = ch_fastqs
             | map { meta, fastqs ->
                   fastqs
             }
             | flatten
             | map { fastq ->
                   def name = fastq.name.toString().replaceAll(/\.trim\.fastq\.gz$/, '') // Remove extension
                   tuple( [id:name, single_end:true], fastq ) // Return tuple
             }
    return ( ch_fastq )
}

// Export the function
return this

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
