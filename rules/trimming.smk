rule trimming:
    """
    Adapter trimming with cutadapt.
    """
    input:
        forward_read = os.path.join(config["raw_reads_dir"], "{sample}_R1_001.fastq.gz"),
        reverse_read = os.path.join(config["raw_reads_dir"], "{sample}_R2_001.fastq.gz")
    output:
        trimmed_forward = "results/trimmed_reads/{sample}/{sample}_R1.fastq.gz",
        trimmed_reverse = "results/trimmed_reads/{sample}/{sample}_R2.fastq.gz"
    params:
        adapter_fwd = config["adapter_fwd"],
        adapter_rev = config["adapter_rev"]
    threads: 4
    log:
        "logs/{sample}/{sample}_trimming.log"
    benchmark:
        "benchmarks/{sample}/{sample}_trimming.txt"
    container:
        "docker://pipecraft/cutadapt:4.4"
    shell:
        """
        cutadapt \
        	-a {params.adapter_fwd} \
        	-A {params.adapter_rev} \
        	-o {output.trimmed_forward} \
        	-p {output.trimmed_reverse} \
			--poly-a \
			--nextseq-trim=20 \
        	--minimum-length=20 \
        	--quality-cutoff=20 \
        	--cores={threads} \
        	{input.forward_read} \
        	{input.reverse_read}
        """