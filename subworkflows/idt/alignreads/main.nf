include { BWAMEM2_MEM }                                 from '../../../modules/idt/bwamem2/mem/main'
include { SENTIEON_BWAMEM }                             from '../../../modules/idt/sentieon/bwamem/main'
include { PICARD_MARKDUPLICATES }                       from '../../../modules/idt/picard/markduplicates/main'
include { SENTIEON_DEDUP }                              from '../../../modules/idt/sentieon/dedup/main'

workflow ALIGNREADS {
    take:
    reads                   // channel: [mandatory] [ meta, reads ]
    fasta                   // channel: [mandatory] [ meta, fasta ]
    fasta_fai               // channel: [mandatory] [ meta, fasta_fai ]
    index                   // path: [mandatory] index/*
    sort_bam                // bool: [mandatory] default: true
    aligner                 // str: "bwa2" or "sentieon-bwamem", default: "bwa2"
    skip_dedup              // bool: [mandatory] default: false

    main:
    // Initialize output channels
    dupmarked_bams = Channel.empty()
    reports = Channel.empty()
    ch_versions = Channel.empty()

    reads_to_align = reads
    switch (aligner) {
        case "bwa2":
            BWAMEM2_MEM(reads,  index, fasta, sort_bam)
            aligned_reads = BWAMEM2_MEM.out.bam
            ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

            if (!skip_dedup) {
                PICARD_MARKDUPLICATES(BWAMEM2_MEM.out.bam, fasta, fasta_fai)
                dupmarked_bams = dupmarked_bams.mix(PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai, failOnDuplicate: true, failOnMismatch: true))
                reports = reports.mix(PICARD_MARKDUPLICATES.out.metrics)
                ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
            }
            break

        case "sentieon-bwamem":
            SENTIEON_BWAMEM(reads, index, fasta, fasta_fai)
            aligned = SENTIEON_BWAMEM.out.bam
            ch_versions = ch_versions.mix(SENTIEON_BWAMEM.out.versions)

            if (!skip_dedup) {
                SENTIEON_DEDUP( SENTIEON_BWAMEM.out.bam_and_bai, fasta, fasta_fai)
                dupmarked_bams = dupmarked_bams.mix(SENTIEON_DEDUP.out.bam.join(SENTIEON_DEDUP.out.bai, failOnDuplicate: true, failOnMismatch: true))
                reports = reports.mix(SENTIEON_DEDUP.out.metrics)
                ch_versions = ch_versions.mix(SENTIEON_DEDUP.out.versions)
            }
            break

        default:
            error "Unknown aligner: ${aligner}"
    }

    emit:
    aligned_bams    = aligned_reads       // channel: [ [meta], bam ]
    dupmarked_bams                        // channel: [ [meta], bam, bai ]
    reports                               // channel: [ [meta], metrics ]

    versions        = ch_versions         // channel: [ ch_versions.yml ]
}
