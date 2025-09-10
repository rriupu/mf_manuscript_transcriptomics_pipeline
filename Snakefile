import os

configfile: "config.yaml"

SAMPLES = glob_wildcards(os.path.join(config["raw_reads_dir"], "{sample}_R1_001.fastq.gz")).sample

include: "rules/qc_before_trimming.smk"
include: "rules/trimming.smk"
include: "rules/qc_after_trimming.smk"
include: "rules/mapping.smk"
# include: "rules/qc_after_mapping.smk"
include: "rules/gene_counting.smk"
# include: "rules/downstream_analysis.smk"

rule all:
    input:
        rules.multiqc_before_trimming.output,
        rules.multiqc_after_trim.output,
        # rules.multiqc_after_alignment.output,
        # expand(rules.htseq_count.output.gene_counts, sample = SAMPLES),
        # rules.data_exploration.output.flag_file,
        # rules.DEA.output.flag_file,
        # rules.gene_clustering.output.flag_file