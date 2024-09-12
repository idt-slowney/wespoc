process YUMI_BAM2ERRORPROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // container "027151828055.dkr.ecr.us-west-2.amazonaws.com/platform-poc-snakemake:latest"

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
    """
    mpileup.py bam2ErrorProfile --inbam $bam \
        --inbed $bed \
        --genome $fasta \
        --thread $task.cpus \
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
