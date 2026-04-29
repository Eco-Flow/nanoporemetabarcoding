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
include { CAT_CAT as CAT_CAT_MEDAKA       } from '../modules/nf-core/cat/cat/main'
//include { DIAMOND_BLASTX                  } from '../modules/nf-core/diamond/blastx/main'
include { BLAST_MAKEBLASTDB               } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN                    } from '../modules/nf-core/blast/blastn/main'
include { SEQKIT_REPLACE                  } from '../modules/nf-core/seqkit/replace/main'
include { BEST_HIT                        } from '../modules/local/blast_best_hit'
include { ASSIGN_TAXONOMY                 } from '../modules/local/assign_taxonomy'
include { paramsSummaryMap                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../subworkflows/local/utils_nfcore_nanoporemetabarcoding_pipeline'
include { validateInputParameters         } from '../subworkflows/local/utils_nfcore_nanoporemetabarcoding_pipeline'

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

    // Prepare the samplesheet channel
    ch_input = ch_samplesheet
             | map {
                meta, fastq ->
                [[id:meta.id, single_end:true], fastq] // Needs to be declared as single end to run PYCHOPPER
             }

    // Prepare metadata channel
    ch_metadata = params.metadata ? Channel.fromPath(params.metadata, checkIfExists: true)
                | splitCsv(header: true)
                | map { row -> [row.id, row.primer_comb, row.sample] }
                | validateInputParameters // Validate metadata so that there are no duplicated values, prob should also check whether fastqs in samplesheet and metadata match
                | map { fastq, primer_comb, sample -> [[id:primer_comb, single_end:true, old_id:fastq], sample] }
                : null // null if no metadata is provided

    // Validation below is working
    ch_metadata_fastqs = ch_metadata.map { meta, sample -> meta.old_id }.distinct().collect().map { ids -> ['validation', ids] }
    ch_input_fastqs    = ch_input.map { meta, fastq -> meta.id }.collect().map { ids -> ['validation', ids] }

    // Validate only if parameter is given
    if (params.metadata) {
        validateSamplesheetMetadata(ch_input_fastqs, ch_metadata_fastqs)
    }

    // Make sure fastqs in samplesheet and metadata match
    //def validateSamplesheetMetadata ( input_channel,metadata_channel )
    //ch_metadata_fastqs.join(ch_input_fastqs)
    //| map { key, input_fastq, metadata_fastq ->
    //                def input_sorted = input_fastq.sort()
    //                def metadata_sorted = metadata_fastq.sort()
    //                        if (metadata_sorted != input_sorted) {
    //                            def missing_in_metadata = input_sorted - metadata_sorted
    //                            def missing_in_input = metadata_sorted - input_sorted
    //
    //                            def error_msg = "ID mismatch between samplesheet and metadata:\n"
    //                            if (missing_in_metadata) error_msg += "In samplesheet but not metadata: ${missing_in_metadata.join(', ')}\n"
    //                            if (missing_in_input) error_msg += "In metadata but not samplesheet: ${missing_in_input.join(', ')}"
    //                            error(error_msg)
    //                        }
    //        return "validation_passed"
    //}
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
        NANOFILT.out.filtreads
    )

    // Flatten the output channel (FASTQs) from cutadapt demultiplex into indidual channels (FASTQ)
    // (check for the function flattenAndMap in the functions section)
    ch_input_f = flattenAndMap(CUTADAPT_F.out.reads, true)

    // Get unkwon reads to reverse complement them later and trim again based on forward barcodes
    // We filter out reads that are unknown as they are probably reverse complemented with regard
    // to the forward barcodes. This is done so that we can trim them again based on the forward barcodes
    // No need for this anymore since cutadapt has an --rc option, so unknown reads are now missing reads
    //ch_unknown = ch_input_f
    //           | filter { meta, fastq ->
    //                    meta.id.contains('unknown')
    //           }

    //
    // MODULE: Run SEQKIT reverse complement
    //

    // Reverse complement reads and run cutadapt again on unknown
    // rc reads and trim again based on forward barcodes
    //SEQKIT_REVCOMP_A (
    //    ch_unknown
    //)

    // Trim based on forward barcodes again
    //CUTADAPT_F_RC (
    //    SEQKIT_REVCOMP_A.out.fastx
    //)

    // Remmove 'unknown_' prefix from metadata and flatten the output channel.
    // This is done so that reads can be concatenated based on the metadata
    // This reads are not unknown anymore
    //ch_unknown = flattenAndMap(CUTADAPT_F_RC.out.reads)
    //           | map { meta, fastq ->
    //               def cleaned_meta = meta.id.replaceFirst(/unknown_/, '') // Remove the 'unknown_' prefix to be able to merge
    //               [meta + [id: cleaned_meta], fastq] // Return updated metadata and fastq
    //           }

    // Group known and unknown (not uknown anymore) reads together based on metadata
    //ch_input_f = ch_input_f
    //           | mix(ch_unknown)
    //           | groupTuple()

    // Concatenate grouped reads together based on metadata.
    //CAT_CAT (
    //    ch_input_f
    //)

    //
    // MODULE: Run SEQKIT reverse complement
    //

    // Barcodes are attached to both ends of the reads, so we need to reverse complement the reads
    // to trim and demultiplex based the other end

    //SEQKIT_REVCOMP_B (
    //    CAT_CAT.out.file_out
    //)

    // Run cutadapt on the reverse complemented reads to trim reverse barcodes
    CUTADAPT_R (
       ch_input_f
    )

    // Filter out FASTQs with less than 10 reads
    ch_input_filtered = CUTADAPT_R.out.reads
                      | flattenAndMap
                      | filter { meta, fastq ->
                           def count = fastq.countFastq()
                           count > params.filt_fastq && !meta.id.contains('unknown') // Filter out FASTQs with less than x reads and with unknown primer combinations
                      }

    //
    // MODULE: Run FastQC
    //
    //FASTQC (
    //    ch_input
    //)
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    //ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    // Prepare data for Nanoplot and amplicon_sorter

    // For running medaka without running minimap2. Medaka already aligns basecalls (amplicons here)
    // to the consensus sequences, so perhaps we can skip minimap2 step, at least for now
    // Make sure this is working as expected

    ch_amplicon_sort = ch_metadata ? ch_input_filtered
                     //| map { meta, fastq -> [meta.id, meta, fastq] } // Extract meta.id replace with sample name according to metadata
                     | join(ch_metadata) // Join on adapter combination
                     | map { meta, fastq, sample ->
                                def new_meta = meta.clone()
                                new_meta.id = sample
                                [new_meta, fastq]
                     }
                     : ch_input_filtered // If no metadata file is provided, use the gunzip output (primer-tag combinations as id)

    //
    // MODULE: Run Gunzip. Might be better to use gunzip from amplicon sorter
    //

    GUNZIP (
        ch_amplicon_sort
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

    // If skip_nanoplot is set to true, skip this module
    ch_nanoplot = params.skip_nanoplot ? Channel.empty() : ch_raw.mix(ch_filt).mix(ch_amplicon_sort)

    NANOPLOT (
        ch_nanoplot
    )

    // // Original name out.txt channel is stats.txt, so multiqc keeps overwritting. Each file needs to have an unique name
    ch_nanoplot_renamed = NANOPLOT.out.txt
                        | map { meta, stats ->
                                    def selectedFile = stats instanceof List ? stats[0] : stats // Not sure why, but sometimes there are two stats files. If that's the case, select the first one [0], which are the original stats
                                    //["${meta.id}_stats.txt", selected_file.text]
                                    def prefix = meta.old_id ? "${meta.old_id}_${meta.id}" :  "${meta.id}"
                                    def renamedFile = selectedFile.copyTo("${workflow.workDir}/renamed_files/${prefix}_stats.txt") //Save renamed files inside the work directory
                                    //def renamed_file = selected_file.copyTo("${meta.id}_stats.txt")
                                    [meta, renamedFile]
                        }

    ch_multiqc_files    = ch_multiqc_files.mix(ch_nanoplot_renamed.map { meta, stats -> stats }).collect()

    //
    // MODULE: Run Amplicon Sorter
    //

    AMPLICON_SORTER (
        GUNZIP.out.gunzip
    )

    // Merge csv from with read counts per ASV. We are doing this to add the read counts to the final table for QC
    ch_merged = AMPLICON_SORTER.out.results_csv
              |collectFile(skip: 2) { meta, csv -> ["${meta.old_id}.merged.txt", csv] } // Merge accoding to old_id (same as the number of tables)
              | map { merged_file ->
                        // Extract the old_id from filename and reconstruct meta
                        def old_id = merged_file.baseName.replaceAll('\\.merged$', '')
                        def new_meta = [id: old_id] // Add other meta fields as needed
                        [new_meta, merged_file]
              }

    //Get group information from amplicon sorter output FASTA files
    ch_group = AMPLICON_SORTER.out.fastas
             | transpose() // Should look this up
             | filter { meta, fasta -> // This is a temporary fix, shoould change this in the module
                            def filename = fasta.name
                            def excludePatterns = ['_unique.fasta', '_consensussequences.fasta', '_nogroup_unique.fasta']
                            return !excludePatterns.any { pattern -> filename.endsWith(pattern) }
             }
             | view()
             | map { meta, fasta ->
                // Get the filename without extension
                 def basename = fasta.baseName
                 // Extract the last part after the final underscore (assuming format from amplicon sorter: name_X_Y)
                 def parts = basename.split('_')
                 def group = "${parts[-2]}_${parts[-1]}" // Get group name from the last two parts of the fasta name without extension

                 def new_meta = meta + [group: group]
                 [new_meta, fasta]
             }

    //
    // MODULE: Run Seqkit Grep
    //

    // Split FASTA files into individual consensus sequence and amplicon sequences
    // for minimap/medaka. Retain group information in the meta so that each consesus
    // sequence is processed together with its respective grouped amplicon sequences
    pattern_amplicons = Channel.of("^\\d+\$") // Amplicon sequences are named with a number
                      | collectFile(name: 'pattern.txt')
    pattern_consensus = Channel.of("^consensus\$") // Consensus sequences are named 'consensus'
                      | collectFile(name: 'pattern.txt')

    // Split into amplicons (each fasta is a group of amplicons without the consensus)
    SEQKIT_AMPLICONS (
        ch_group,
        pattern_amplicons.first()
    )

    GUNZIP_SEQKIT_GREP_A (
        SEQKIT_AMPLICONS.out.filter
    )

    // Split into consensus (each fasta is a consensus without the grouped amplicons)
    SEQKIT_CONSENSUS (
        ch_group,
        pattern_consensus.first()
    )

    // Rename the consensus sequences to their group and meta id
    SEQKIT_REPLACE(
        SEQKIT_CONSENSUS.out.filter
    )

    GUNZIP_SEQKIT_GREP_C (
        SEQKIT_REPLACE.out.fastx
    )

    // Join consensus and amplicon sequences based on metadata and separate them in
    // a multichannel (keeps grouped amplicons and their respective consensus sequence in sync)
    //ch_minimap = GUNZIP_SEQKIT_GREP_A.out.gunzip // Grouped amplicons
    //           | join(GUNZIP_SEQKIT_GREP_C.out.gunzip) // Consensus
    //           | multiMap { meta, amps, cons -> // meta: metadata, amps: amplicon sequences, cons: consensus sequences
    //                        amps : [ meta, amps] // Return a tuple with metadata and amplicon sequences
    //                        cons : [ meta, cons ] // Return a tuple with metadata and consensus sequences
    //           }

    // For running medaka without running minimap2. Medaka already aligns basecalls (amplicons here)
    // to the consensus sequences, so perhaps we can skip minimap2 step, at least for now
    // Make sure this is working as expected
    ch_medaka = GUNZIP_SEQKIT_GREP_A.out.gunzip
              | join(GUNZIP_SEQKIT_GREP_C.out.gunzip)

    //
    // MODULE: Run Minimap2
    //



    //
    // MODULE: Run Medaka
    //


    MEDAKA (
        ch_medaka
    )

    // Concatenate corrected consensus sequences so they can be blasted all together
    // This is very important buecause if they are blasted induvidually the dabaase has
    // to be loaded into memory every time
    ch_corrected = MEDAKA.out.assembly
                 | map {
                    meta, fasta ->
                    [[id:meta.old_id], fasta] // old_id to concatenate them
                 }
                 | groupTuple(by: 0)

    CAT_CAT_MEDAKA (
        ch_corrected
    )

    ch_corrected_concat = CAT_CAT_MEDAKA.out.file_out

    //
    // MODULE: Run makeblastdb
    //

    // Prepare ch_databse channel to build a custom database for blast
    ch_database = (params.custom_db && !params.blast_db) ? Channel.fromPath(params.custom_db)
                | map { db ->
                        [[id:'database'], db]
                } : Channel.empty() // If no custom database is provided, use an empty channel


    BLAST_MAKEBLASTDB (
        ch_database
    )

    // Mix in case an already built blast database is already give. People should only input a
    // path to make the database or give the built database. This shouldn't be possible,
    // will write code later to prevent it
    ch_prebuilt_db = params.blast_db ? Channel.fromPath(params.blast_db)
                   | map { db ->
                            [[id:'prebuilt_database'], db]
                   } : Channel.empty() // If no prebuild database is provided, use an empty channel


    ch_blast = BLAST_MAKEBLASTDB.out.db.mix(ch_prebuilt_db)

    //
    // MODULE: Run BLAST
    //

    BLAST_BLASTN (
        ch_corrected_concat,
        ch_blast.first() // .first() so that channel can be used several times
    )

    //
    // MODULE: Run blast best hit
    //

    // To annotate the consensus (or ASVs) sequences, we need to estabish a criteria
    // to select the best blast hit. In this case the best blast hit is established first
    // by bitscore and second by evalue

    //BEST_HIT (
    //    BLAST_BLASTN.out.txt
    //)

    // Join blast hits channel with number of reads per ASV channel
    ch_assign_taxonomy = BLAST_BLASTN.out.txt
                       | join(ch_merged, by: 0)
                       | multiMap { meta, blast, csv ->
                                blast: [meta, blast]
                                csv: [meta,csv]
                       }

    ch_assign_taxonomy.blast.view()
    ch_assign_taxonomy.csv.view()

    //
    // MODULE: Run assign taxonomy
    //
    ch_sql_db = Channel.fromPath(params.sql_db)

    ASSIGN_TAXONOMY(
        ch_assign_taxonomy.blast,
        ch_assign_taxonomy.csv,
        ch_sql_db.first()
    )

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

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
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

    // Make sure fastqs in samplesheet and metadata match
    def validateSamplesheetMetadata (input_channel, metadata_channel) {
        metadata_channel
            |join(input_channel)
            | map { key, input_fastq, metadata_fastq ->
                        def input_sorted = input_fastq.sort()
                        def metadata_sorted = metadata_fastq.sort()
                                if (metadata_sorted != input_sorted) {
                                    def missing_in_metadata = input_sorted - metadata_sorted
                                    def missing_in_input = metadata_sorted - input_sorted

                                    def error_msg = "ID mismatch between samplesheet and metadata:\n"
                                    if (missing_in_metadata) error_msg += "In metadata but not samplesheet: ${missing_in_metadata.join(', ')}\n"
                                    if (missing_in_input) error_msg += "In samplesheet but not metadata: ${missing_in_input.join(', ')}"
                                    error(error_msg)
                                }
                        "validation_passed"
            }
    }

def flattenAndMap(ch_fastqs, preserve_old_id = false) { // preserve_old_id becuase this change only needs to be done on the first cutadapt process
    ch_fastq = ch_fastqs
             | flatMap { meta, fastqs ->
                   // Use flatMap instead of map + flatten
                   fastqs.collect { fastq ->
                       def name = fastq.name.toString().replaceAll(/\.trim\.fastq\.gz$/, '')
                       def new_meta = meta + [id: name]
                       // Only add old_id if requested
                       if (preserve_old_id) {
                           new_meta = new_meta + [old_id: meta.id]
                       }
                       tuple(new_meta, fastq)
                   }
             }

    return ch_fastq
}

// Export the function
return this

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
