/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                          } from '../modules/nf-core/fastqc/main'
include { MULTIQC                         } from '../modules/nf-core/multiqc/main'
include { PORECHOP_PORECHOP               } from '../modules/nf-core/porechop/porechop/main'
include { NANOPLOT                        } from '../modules/nf-core/nanoplot/main'
include { NANOFILT                        } from '../modules/nf-core/nanofilt/main'
include { SEQKIT_SEQ as SEQKIT_REVCOMP_A  } from '../modules/nf-core/seqkit/seq/main'
include { SEQKIT_SEQ as SEQKIT_REVCOMP_B  } from '../modules/nf-core/seqkit/seq/main'
include { CUTADAPT as CUTADAPT_F          } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_F_RC       } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_R          } from '../modules/nf-core/cutadapt/main'
include { AMPLICON_SORTER                 } from '../modules/local/amplicon_sorter'
include { SEQKIT_GREP as SEQKIT_AMPLICONS } from '../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_CONSENSUS } from '../modules/nf-core/seqkit/grep/main'
include { GUNZIP                          } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_SEQKIT_GREP_A  } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_SEQKIT_GREP_C  } from '../modules/nf-core/gunzip/main'
include { MEDAKA                          } from '../modules/nf-core/medaka/main'
include { CAT_CAT                         } from '../modules/nf-core/cat/cat/main'
include { paramsSummaryMap                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../subworkflows/local/utils_nfcore_nanoporemetabarcoding_pipeline'

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

    // Flatten the output channel (FASTQs) from cutadapt demultiplex into indidual channels (FASTQ)
    // (check for the function flattenAndMap in the functions section)
    ch_input_f = flattenAndMap(CUTADAPT_F.out.reads)

    // Get unkwon reads to reverse complement them later and trim again based on forward barcodes
    // We filter out reads that are unknown as they are probably reverse complemented with regard
    // to the forward barcodes. This is done so that we can trim them again based on the forward barcodes
    ch_unknown = ch_input_f
               | filter { meta, fastq ->
                        meta.id.contains('unknown')
               }

    //
    // MODULE: Run SEQKIT reverse complement
    //

    // Reverse complement reads and run cutadapt again on unknown
    // rc reads and trim again based on forward barcodes
    SEQKIT_REVCOMP_A (
        ch_unknown
    )

    // Trim based on forward barcodes again
    CUTADAPT_F_RC (
        SEQKIT_REVCOMP_A.out.fastx
    )

    // Remmove 'unknown_' prefix from metadata and flatten the output channel.
    // This is done so that reads can be concatenated based on the metadata
    ch_unknown = flattenAndMap(CUTADAPT_F_RC.out.reads)
               | map { meta, fastq ->
                   def cleaned_meta = meta.id.replaceFirst(/unknown_/, '') // Remove the 'unknown_' prefix to be able to merge
                   [[id: cleaned_meta, single_end: meta.single_end], fastq] // Return updated metadata and fastq
               }

    // Group known and unknown (not uknown anymore) reads together based on metadata
    ch_input_f = ch_input_f
               | mix(ch_unknown)
               | groupTuple()

    // Concatenate grouped reads together based on metadata.
    CAT_CAT (
        ch_input_f
    )

    //
    // MODULE: Run SEQKIT reverse complement
    //

    // Barcodes are attached to both ends of the reads, so we need to reverse complement the reads
    // to trim and demultiplex based the other end

    SEQKIT_REVCOMP_B (
        CAT_CAT.out.file_out
    )

    // Run cutadapt on the reverse complemented reads to trim reverse barcodes
    CUTADAPT_R (
       SEQKIT_REVCOMP_B.out.fastx 
    )

    // Filter out FASTQs with less than 10 reads
    ch_input_filtered = CUTADAPT_R.out.reads
                      | flattenAndMap
                      | filter { meta, fastq ->
                           def count = fastq.countFastq()
                           count > 100 && !meta.id.contains('unknown') // Filter out FASTQs with less than 1000 reads
                      }

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
        ch_raw.mix(ch_filt).mix(ch_input_filtered)
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
    // MODULE: Run Gunzip. Might be better to use gunzip from amplicon sorter
    //

    GUNZIP (
        ch_input_filtered
    )
 
    //
    // MODULE: Run Amplicon Sorter
    //

    AMPLICON_SORTER (
        GUNZIP.out.gunzip
    )
    
    //Get group information from amplicon sorter output FASTA files
    ch_group = AMPLICON_SORTER.out.fastas
             | transpose()
             | map { meta, fasta ->
                // Get the filename without extension
                 def basename = fasta.baseName
                 // Extract the last part after the final underscore (assuming format from amplicon sorter: name_X_Y)
                 def parts = basename.split('_')
                 def group = "${parts[-2]}_${parts[-1]}" // Get group name from the last two parts of the fasta name without extension

                 def new_meta = meta + [group: group]
                 [new_meta, fasta]
             }
    
    ch_group.view()

    //
    // MODULE: Run Seqkit Grep
    //

    // Split FASTA files into individual into consensus sequences and amplicon sequences
    // for minimap/medaka. Retain group information in the meta
    pattern_amplicons = Channel.of("^\\d+\$") // Amplicon sequences are named with a number
                      | collectFile(name: 'pattern.txt')
    pattern_consensus = Channel.of("^consensus\$") // Consensus sequences are named 'consensus'
                      | collectFile(name: 'pattern.txt')

    SEQKIT_AMPLICONS (
        ch_group,
        pattern_amplicons.first()
    )

    GUNZIP_SEQKIT_GREP_A (
        SEQKIT_AMPLICONS.out.filter
    )

    SEQKIT_CONSENSUS (
        ch_group,
        pattern_consensus.first()
    )

    GUNZIP_SEQKIT_GREP_C (
        SEQKIT_CONSENSUS.out.filter
    )

    // Join consensus and amplicon sequences based on metadata and separate them in
    // a multichannel (keeps them in sync but can be processed separately)
    ch_minimap = GUNZIP_SEQKIT_GREP_A.out.gunzip
               | join(GUNZIP_SEQKIT_GREP_C.out.gunzip)
               | multiMap { meta, amps, cons -> // meta: metadata, amps: amplicon sequences, cons: consensus sequences
                            amps : [ meta, amps] // Return a tuple with metadata and amplicon sequences
                            cons : [ meta, cons ] // Return a tuple with metadata and consensus sequences
                }
    // For running medaka without running minimap2. Medaka already aligns basecalls (amplicons here)
    // to the consensus sequences, so perhaps we can skip minimap2 step
    ch_medaka = GUNZIP_SEQKIT_GREP_A.out.gunzip
              | join(GUNZIP_SEQKIT_GREP_C.out.gunzip)

    //
    // MODULE: Run Minimap2
    //

    MEDAKA (
        ch_medaka
    )

    //
    // MODULE: Run Medaka
    //



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
