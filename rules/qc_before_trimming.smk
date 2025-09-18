rule fastqc_before_trimming:
    """
    Run FastQC on raw reads.
    """
    input:
        forward_read = os.path.join(config["raw_reads_dir"], "{sample}_R1_001.fastq.gz"),
        reverse_read = os.path.join(config["raw_reads_dir"], "{sample}_R2_001.fastq.gz")
    output:
        html1 = os.path.join(config["output_dir"], "qc_reports", "before_trimming", "{sample}", "{sample}_R1_001_fastqc.html"),
        zipfile1 = os.path.join(config["output_dir"], "qc_reports", "before_trimming", "{sample}", "{sample}_R1_001_fastqc.zip"),
        html2 = os.path.join(config["output_dir"], "qc_reports", "before_trimming", "{sample}", "{sample}_R2_001_fastqc.html"),
        zipfile2 = os.path.join(config["output_dir"], "qc_reports", "before_trimming", "{sample}", "{sample}_R2_001_fastqc.zip")
    params:
        wd = os.path.join(config["output_dir"], "qc_reports", "before_trimming", "{sample}/")
    threads: 4
    log:
        os.path.join(config["logs_dir"], "FastQC", "before_trimming", "{sample}", "{sample}_fastqc_before_trimming.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "FastQC", "before_trimming", "{sample}", "{sample}_fastqc_before_trimming.txt")
    container:
        "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    shell:
        """
        mkdir -p {params.wd} 2>> {log}
        
        fastqc \
            --format fastq \
            --threads {threads} \
            --outdir {params.wd} \
            --dir {params.wd} \
            {input.forward_read} {input.reverse_read} \
        &>> {log}
        """


rule multiqc_before_trimming:
    """
    Run mulqiQC on fastqc reports.
    """
    input:
        expand(rules.fastqc_before_trimming.output.html1, sample = SAMPLES)
    output:
        report = os.path.join(config["output_dir"], "qc_reports", "before_trimming", "multiqc_report.html")
    params:
        wd = os.path.join(config["output_dir"], "qc_reports", "before_trimming/")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "multiqc", "before_trimming.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "multiqc", "before_trimming.txt")
    container:
        "docker://multiqc/multiqc:v1.29"
    shell:
        """
        cd {params.wd}

        multiqc .
        """