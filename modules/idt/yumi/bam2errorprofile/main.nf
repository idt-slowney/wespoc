process YUMI_BAM2ERRORPROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "slowney/yumi:latest"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bed)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*_mpileup_error_rate.txt"), emit: error_rate
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tmpDir = params.tmpDir ?: "/tmp"
    """
    mpileup.py bam2ErrorProfileParallel --inbam $bam \
        --inbed $bed \
        --genome $fasta \
        --thread $task.cpus \
        --tmpDir $tmpDir \
        $args \
        --outErrorProfile ${prefix}_mpileup_error_rate.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam2ErrorProfile: \$(python --version |& sed '1!d ; s/Python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mpileup_error_rate.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mpileup: \$(python --version |& sed '1!d ; s/Python //')
    END_VERSIONS
    """
}
