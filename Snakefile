import os

configfile: "config.yaml"

SAMPLES = glob_wildcards(os.path.join(config["raw_reads_dir"], "{sample}_R1_001.fastq.gz")).sample

include: "rules/qc_before_trimming.smk"
include: "rules/trimming.smk"
include: "rules/qc_after_trimming.smk"
include: "rules/mapping.smk"
include: "rules/gene_counting.smk"
include: "rules/downstream_analysis.smk"
include: "rules/unibind_TFBS_enrichment.smk"

rule all:
    input:
        rules.multiqc_before_trimming.output,
        rules.multiqc_after_trim.output,
        rules.data_exploration.output.sample_sample_distance_heatmap,
        # rules.gene_clustering.output.clusters,
        # rules.GSVA.output.GSVA_plot,
        expand(rules.generate_TFBS_enrichment_summary_plots.output.summary_plot, deg_subset = ["all", "upregulated", "downregulated"])
        # rules.DEA.output.flag_file,