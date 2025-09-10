rule get_genome:
    """
    Get the genome fasta file.
    """
    output:
        genome_fasta = os.path.join(config["data_dir"], "genome", os.path.splitext(os.path.basename(config["genome_url"]))[0])
    params:
        genome_url = config["genome_url"]
    threads: 1
    log:
        os.path.join(config["logs_dir"], "genome_download", "download.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "genome_download", "download.txt")
    shell:
        """
        wget -O - {params.genome_url} \
        | gunzip -c > {output.genome_fasta} \
        2> {log}
        """


rule get_gtf:
    """
    Get GTF file with gene annotations.
    """
    output:
        gtf_file = os.path.join(config["data_dir"], "annotations", "hg38.ncbiRefSeq.gtf")
    params:
        gtf_url = config["gtf_url"]
    threads: 1
    log:
        os.path.join(config["logs_dir"], "get_gtf", "get_gtf.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "get_gtf", "get_gtf.txt")
    shell:
        """
        wget -O - {params.gtf_url} \
        | gunzip -c > {output.gtf_file} \
        2> {log}
        """


rule genome_indexing:
    """
    Index the genome using STAR.
    """
    input:
        genome_fasta = rules.get_genome.output.genome_fasta,
        gtf_file = rules.get_gtf.output.gtf_file
    output:
        touch("data/genome/indexing_done")
    params:
        genome_dir = os.path.join(config["data_dir"], "star_index")
    threads: 10
    log:
        os.path.join(config["logs_dir"], "genome_indexing", "indexing.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "genome_indexing", "indexing.log")
    container:
        "docker://mgibio/star:2.7.0f"
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {params.genome_dir} \
            --genomeFastaFiles {input.genome_fasta} \
            --sjdbGTFfile {input.gtf_file} \
            --sjdbOverhang ReadLength-1 \
        2> {log}
        """

rule star_alignment:
    """
    Align reads with STAR.
    """
    input:
        trimmed_forward = rules.trimming.output.trimmed_forward,
        trimmed_reverse = rules.trimming.output.trimmed_reverse,
        genome_indexing_flag_file = rules.genome_indexing.output
    output:
        aligned_reads = os.path.join("results", "aligned_reads", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
    params:
        genome_dir = os.path.join(config["data_dir"], "star_index")
    threads: 8
    log:
        os.path.join(config["logs_dir"], "mapping", "{sample}", "{sample}_star_mapping.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "mapping", "{sample}", "{sample}_star_mapping.log")
    container:
        "docker://mgibio/star:2.7.0f"
    shell:
        """
        STAR \
			--readFilesCommand gunzip -c \
			--genomeDir {params.genome_dir} \
			--runThreadN {threads} \
			--readFilesIn {input.trimmed_forward} {input.trimmed_reverse} \
			--outFileNamePrefix $(dirname {output.aligned_reads}) \
			--outSAMtype BAM SortedByCoordinate \
			--outSAMunmapped Within \
			--outSAMattributes Standard \
        2> {log}
        """
