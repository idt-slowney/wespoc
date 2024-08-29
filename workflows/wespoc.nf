/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_wespoc_pipeline'

include { MULTIQC                } from '../modules/nf-core/multiqc/main'

// Custom Subworkflows
include { PROCESSFASTQ              } from '../subworkflows/idt/processfastq/main'
include { ALIGNREADS                } from '../subworkflows/idt/alignreads/main'
include { COLLECT_QC_METRICS        } from '../subworkflows/idt/collectqcmetrics/main'

// Custom Modules
include { VARDICTJAVA               } from '../modules/idt/vardictjava/main'
include { YUMI_BAM2ERRORPROFILE     } from '../modules/idt/yumi/bam2errorprofile/main'
include { PICARD_BEDTOINTERVALLIST  } from '../modules/idt/picard/bedtointervallist/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow wespoc {
    take:
        reads
        fasta
        fasta_dict
        fasta_fai
        index 
        aligner
        target_bed

    main:
    // Instantiate empty channels
    ch_multiqc_files = Channel.empty()
    versions = Channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        YOUR WORKFLOW HERE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // Prep Intervals 
    PICARD_BEDTOINTERVALLIST( target_bed,
                              fasta_dict,
                              [] )

    // Trimming, downsampling, FastQC
    PROCESSFASTQ ( reads, 
                   params.skip_fastqc, 
                   params.skip_downsample, 
                   params.umi_flag, 
                   params.skip_extractumi, 
                   params.read_structure,
                   params.skip_trimming )
    versions = versions.mix(PROCESSFASTQ.out.versions)

    // Align and dedup
    ALIGNREADS ( PROCESSFASTQ.out.reads, 
                 fasta, 
                 fasta_fai, 
                 index, 
                 params.sort_bam, 
                 aligner, 
                 params.skip_dedup )
    versions = versions.mix(ALIGNREADS.out.versions)

    // QC and Error profiling
    interval_list = PICARD_BEDTOINTERVALLIST.out.interval_list
                        .map{ meta, interval_list -> interval_list }
                        // .flatten()

    COLLECT_QC_METRICS( ALIGNREADS.out.dupmarked_bams, 
                        interval_list,
                        interval_list,
                        fasta,
                        fasta_fai,
                        fasta_dict )
    versions = versions.mix(COLLECT_QC_METRICS.out.versions)

    YUMI_BAM2ERRORPROFILE( ALIGNREADS.out.dupmarked_bams.map{meta, bam, bai -> [meta, bam]}, 
                           target_bed, 
                           fasta )
    versions = versions.mix(YUMI_BAM2ERRORPROFILE.out.versions)

    // Variant Calling
    varCallIn_ch = ALIGNREADS.out.dupmarked_bams
                        .combine(target_bed.map{meta, bed -> bed})

    VARDICTJAVA(varCallIn_ch, fasta, fasta_fai)

    
    
    // versions = versions.mix(VARDICTJAVA.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF YOUR WORKFLOW
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
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

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
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
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
