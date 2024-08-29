#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IDT/wespoc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/IDT/wespoc
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { wespoc  } from './workflows/wespoc'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_wespoc_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_wespoc_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_wespoc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta                   = getGenomeAttribute('fasta')
params.fasta_dict              = getGenomeAttribute('fasta_dict')
params.fasta_fai               = getGenomeAttribute('fasta_fai')
params.bwa                     = getGenomeAttribute('bwa')
params.bwamem2                 = getGenomeAttribute('bwamem2')
// params.target_bed              = Channel.fromPath("/mnt/archive/work/slowney/Genomes/beds/Exome_Exomev2_Targets.bed").collect()

include { PREPREFERENCEGENOME     } from './subworkflows/idt/prepreferencegenome/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALIZE PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

aligner                 = params.aligner
fasta                   = params.fasta              ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() 
                                                    : Channel.empty()
fasta_fai               = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect()            
                                                    : Channel.empty()
fasta_dict              = params.fasta_dict         ? Channel.fromPath(params.fasta_dict).collect()           
                                                    : Channel.empty()
target_bed              = params.target_bed         ? Channel.fromPath(params.target_bed).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                    : Channel.value([[],[]])
//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow IDT_wespoc {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    versions = Channel.empty()

    PREPREFERENCEGENOME(fasta)

    // Gather any new genome files from PREPREFERENCEGENOME that aren't specified in params
    bwa         = params.bwa        ? Channel.fromPath(params.bwa).map{ it -> [ [id:'bwa'], it ] }.collect()
                                    : PREPREFERENCEGENOME.out.bwa
    bwamem2     = params.bwamem2    ? Channel.fromPath(params.bwamem2).map{ it -> [ [id:'bwamem2'], it ] }.collect()
                                    : PREPREFERENCEGENOME.out.bwamem2
    fasta_dict  = params.fasta_dict ? Channel.fromPath(params.fasta_dict).map{ it -> [ [id:'fasta_dict'], it ] }.collect()
                                    : PREPREFERENCEGENOME.out.fasta_dict
    fasta_fai   = params.fasta_fai  ? Channel.fromPath(params.fasta_fai).map{ it -> [ [id:'fai'], it ] }.collect()
                                    : PREPREFERENCEGENOME.out.fasta_fai
    // Generalize index var
    index = (aligner == "sentieon-bwamem") ? bwa : bwamem2


    versions = versions.mix(PREPREFERENCEGENOME.out.versions)

    //
    // WORKFLOW: Run pipeline
    //
    wespoc (
        samplesheet,
        fasta,
        fasta_dict,
        fasta_fai,
        index,
        aligner,
        target_bed
    )

    emit:
    multiqc_report = wespoc.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    IDT_wespoc (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        IDT_wespoc.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
