rule fastqc_after_alignment:
    """
    Run FastQC on all samples.
    """
    input:
        aligned_reads = rules.star_alignment.output.aligned_reads
    output:
        html = "results/qc_reports/after_alignment/{sample}/aligned_reads_fastqc.html",
        zipfile = "results/qc_reports/after_alignment/{sample}/aligned_reads_fastqc.zip"
    params:
        wd = "results/qc_reports/after_alignment/{sample}/"
    threads: 4
    log:
        "logs/{sample}/{sample}_fastqc_after_alignment.log"
    benchmark:
        "benchmarks/{sample}/{sample}_fastqc_after_alignment.txt"
    container:
        "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    shell:
        """
        mkdir -p {params.wd} 2>> {log}
        
        fastqc \
            --format bam \
            --threads {threads} \
            --outdir {params.wd} \
            --dir {params.wd} \
            {input.aligned_reads} \
        &>> {log}
        """

rule multiqc_after_alignment:
    """
    Run mulqiQC on fastqc reports.
    """
    input:
        expand(rules.fastqc_after_alignment.output.html, sample = SAMPLES)
    output:
        report = "results/qc_reports/after_alignment/multiqc_report.html"
    params:
        wd = "results/qc_reports/after_alignment/"
    threads: 1
    log:
        "logs/multiqc/after_alignment.log"
    benchmark:
        "benchmarks/multiqc/after_alignment.txt"
    container:
        "docker://multiqc/multiqc:v1.29"
    shell:
        """
        cd {params.wd}

        multiqc .
        """