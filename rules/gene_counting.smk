rule htseq_count:
    """
    Get gene counts using HTSeq.
    """
    input:
        gtf_file = rules.get_gtf.output.gtf_file,
        aligned_reads = rules.star_alignment.output.aligned_reads
    output:
        gene_counts = os.path.join("results", "gene_counts", "{sample}_geneCounts.tsv")
    log:
        "logs/{sample}/{sample}_htseq_count.log"
    benchmark:
        "benchmarks/{sample}/{sample}_htseq_count.txt"
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
            {input.gtf_file}
        """
