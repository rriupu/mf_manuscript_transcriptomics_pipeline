rule htseq_count:
    """
    Get gene counts using HTSeq.
    """
    input:
        gtf_file = rules.get_gtf.output.gtf_file,
        aligned_reads = rules.star_alignment.output.aligned_reads
    output:
        gene_counts = os.path.join(config["output_dir"], "gene_counts", "{sample}_geneCounts.tsv")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "htseq", "{sample}", "{sample}_htseq_count.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "htseq", "{sample}", "{sample}_htseq_count.txt")
    container:
        "docker://pegi3s/htseq:2.0.9"
    shell:
        """
        htseq-count \
            -s yes \
            -r pos \
            -t exon \
            -i gene_id \
            -c {output.gene_counts} \
            {input.aligned_reads} \
            {input.gtf_file} \
        2> {log}
        """
