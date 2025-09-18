rule fastqc_after_trim:
    """
    Run FastQC on all samples.
    """
    input:
        trimmed_forward = rules.trimming.output.trimmed_forward,
        trimmed_reverse = rules.trimming.output.trimmed_reverse
    output:
        html1 = os.path.join(config["output_dir"], "qc_reports", "after_trimming", "{sample}", "{sample}_R1_fastqc.html"),
        zipfile1 = os.path.join(config["output_dir"], "qc_reports", "after_trimming", "{sample}", "{sample}_R1_fastqc.zip"),
        html2 = os.path.join(config["output_dir"], "qc_reports", "after_trimming", "{sample}", "{sample}_R2_fastqc.html"),
        zipfile2 = os.path.join(config["output_dir"], "qc_reports", "after_trimming", "{sample}", "{sample}_R2_fastqc.zip")
    params:
        wd = os.path.join(config["output_dir"], "qc_reports", "after_trimming", "{sample}/")
    threads: 5
    log:
        os.path.join(config["logs_dir"], "FastQC", "after_trimming", "{sample}", "{sample}_fastqc_after_trimming.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "FastQC", "after_trimming", "{sample}", "{sample}_fastqc_after_trimming.txt")
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
            {input.trimmed_forward} \
            {input.trimmed_reverse} \
        &>> {log}
        """

rule multiqc_after_trim:
    """
    Run mulqiQC on fastqc reports.
    """
    input:
        expand(rules.fastqc_after_trim.output.html1, sample = SAMPLES)
    output:
        report = os.path.join(config["output_dir"], "qc_reports", "after_trimming", "multiqc_report.html")
    params:
        wd = os.path.join(config["output_dir"], "qc_reports", "after_trimming/")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "multiqc", "after_trimming.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "multiqc", "after_trimming.txt")
    container:
        "docker://multiqc/multiqc:v1.29"
    shell:
        """
        cd {params.wd}

        multiqc .
        """