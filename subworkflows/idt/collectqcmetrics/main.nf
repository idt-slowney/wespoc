include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/idt/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS       } from '../../../modules/idt/picard/collecthsmetrics/main'
include { SENTIEON_DATAMETRICS          } from '../../../modules/idt/sentieon/datametrics/main'
// include { SENTIEON_HSMETRICS            } from '../../../modules/idt/sentieon/hsmetrics/main'
include { PICARD_COLLECTOXOGMETRICS     } from '../../../modules/idt/picard/collectoxogmetrics/main'
include { PICARD_COLLECTSEQUENCINGARTIFACTMETRICS } from '../../../modules/idt/picard/collectsequencingartifactmetrics/main'

workflow COLLECT_QC_METRICS {
    take:
    bam_and_bai  // channel: [ val(meta), [bam], [bai] ]
    bait_interval
    target_interval
    fasta                // channel: [ val(meta), fasta ]
    fasta_fai            // channel: [ val(meta), fasta_fai ]
    fasta_dict           // channel: [ val(meta), fasta_dict ]

    main:
    versions = Channel.empty()
    coverage_metrics = Channel.empty()
    multiple_metrics = Channel.empty()
    oxog_metrics = Channel.empty()
    artifact_metrics = Channel.empty()

    // -------------------------------- MultipleMetrics -------------------------------- //
    // GCBias, MeanQualityByCycle, QualityDistribution, InsertSizeMetrics, AlignmentStats
    PICARD_COLLECTMULTIPLEMETRICS( bam_and_bai, fasta, fasta_fai ) // ext.when = { accelerate_metrics_flag == false }
    // SENTIEON_DATAMETRICS( bam_and_bai, fasta, fasta_fai ) // ext.when = { accelerate_metrics_flag == true }

    multiple_metrics = multiple_metrics.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics)
    // multiple_metrics = multiple_metrics.mix(SENTIEON_DATAMETRICS.out.metrics)
    versions = versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    // versions = versions.mix(SENTIEON_DATAMETRICS.out.versions.first())

    // ----------------------------------- HsMetrics ----------------------------------- //
    // ch_hsMetrics = bam_and_bai.map{ meta, bam, bai -> [ meta, bam, bai, bait_interval, target_interval ] }
    ch_hsMetrics = bam_and_bai.combine(bait_interval).combine(target_interval)
    ch_hsMetrics.view()

    PICARD_COLLECTHSMETRICS( ch_hsMetrics, fasta, fasta_fai, fasta_dict )
    // SENTIEON_HSMETRICS()

    versions = versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
    // versions = versions.mix(SENTIEON_HSMETRICS.out.versions.first())
    coverage_metrics = coverage_metrics.mix(PICARD_COLLECTHSMETRICS.out.metrics)
    // coverage_metrics = coverage_metrics.mix(SENTIEON_HSMETRICS.out.metrics)

    // ----------------------------------- SequenceArtifacts ----------------------------------- //
    PICARD_COLLECTOXOGMETRICS( bam_and_bai.map{ meta, bam, bai -> [meta, bam] }, fasta)
    oxog_metrics = oxog_metrics.mix(PICARD_COLLECTOXOGMETRICS.out.metrics)
    versions = versions.mix(PICARD_COLLECTOXOGMETRICS.out.versions)

    // ----------------------------------- OxoMetrics ----------------------------------- //
    PICARD_COLLECTSEQUENCINGARTIFACTMETRICS( bam_and_bai.map{ meta, bam, bai -> [meta, bam] }, fasta )
    artifact_metrics = artifact_metrics.mix( PICARD_COLLECTSEQUENCINGARTIFACTMETRICS.out.metrics )
    versions = versions.mix( PICARD_COLLECTSEQUENCINGARTIFACTMETRICS.out.versions )

    emit:
    coverage_metrics    = coverage_metrics                       // channel: [ val(meta), [ coverage_metrics ] ]
    multiple_metrics    = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), [ multiple_metrics ] ]
    oxog_metrics        = oxog_metrics
    artifact_metrics    = artifact_metrics
    versions            = versions                               // channel: [ versions.yml ]
}
