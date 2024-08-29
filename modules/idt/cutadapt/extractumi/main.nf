process CUTADAPT_EXTRACTUMI {
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads)
    val(read_structure)

    output:
    tuple val(meta), path("*.ExtractedUMI.fastq.gz"), emit: reads
    path  "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cores = task.cpus ?: 1
    def prefix = task.ext.prefix ?: ""
    def ReadStructure = read_structure.count("M") == 1 ? read_structure.tokenize("M")[0] : [read_structure.tokenize(" ")[0].tokenize("M")[0],read_structure.tokenize(" ")[1].tokenize("M")[0]]
    if (meta.single_end) {
        def args_list = args.split("\\s(?=--)").toList()
        args_list.removeAll { it.toLowerCase().contains('_r2 ') }
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        cutadapt \\
            $args \\
            -u ${ReadStructure} \\
            --rename='{id}\\tZA:Z:{r1.cut_prefix}' \\
            --cores $cores \\
            -o ${prefix}.ExtractedUMI.fastq.gz \\
            ${reads[0]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_R1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_R1.fastq.gz
        [ ! -f  ${prefix}_R2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_R2.fastq.gz
        cutadapt \\
            $args \\
            -u ${ReadStructure[0]} -U ${ReadStructure[1]} \\
            --rename='{id}\\tZA:Z:{r1.cut_prefix}\\tZB:Z:{r2.cut_prefix}\\tRX:Z:{r1.cut_prefix}-{r2.cut_prefix}' \\
            --cores $cores \\
            -o ${prefix}_R1.ExtractedUMI.fastq.gz -p ${prefix}_R2.ExtractedUMI.fastq.gz \\
            ${reads[0]} ${reads[1]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }
}
