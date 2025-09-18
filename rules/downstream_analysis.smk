rule data_exploration:
    """
    Run an initial data exploratory analysis with R.
    """
    input:
        gene_counts = expand(rules.htseq_count.output.gene_counts, sample = SAMPLES),
        condition_mapping = "conditionMapping.tsv"
    output:
        sample_sample_distance_heatmap = os.path.join(config["output_dir"], "downstream_analysis", "data_exploration", "sample_sample_distance_blind_no_donor_effect.pdf"),
        correlation_heatmap = os.path.join(config["output_dir"], "downstream_analysis", "data_exploration", "correlation_heatmap_no_donor_effect.pdf"),
        PCA_plot = os.path.join(config["output_dir"], "downstream_analysis", "data_exploration", "PCA_blind_no_donor_effect.pdf"),
        PCA_plot_facet = os.path.join(config["output_dir"], "downstream_analysis", "data_exploration", "PCA_blind_no_donor_effect_facet.pdf")
    params:
        r_script = os.path.join(config["scripts_dir"], "data_exploration.R"),
        gene_counts_dir = os.path.join(config["output_dir"], "gene_counts"),
        output_dir = os.path.join(config["output_dir"], "downstream_analysis", "data_exploration")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "data_exploration", "data_exploration.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "data_exploration", "data_exploration.txt")
    conda:
        "../envs/downstream.yaml"
    shell:
        """
        Rscript {params.r_script} \
            --conditionMapping {input.condition_mapping} \
            --geneCountsDir {params.gene_counts_dir} \
            --outputDir {params.output_dir}
        """

checkpoint DEA:
    """
    Run a differential expression analysis with R.
    """
    input:
        gene_counts = expand(rules.htseq_count.output.gene_counts, sample = SAMPLES),
        condition_mapping = "conditionMapping.tsv",
        r_script = os.path.join(config["scripts_dir"], "DEA.R")
    output:
        comparison_results = directory(os.path.join(config["output_dir"], "downstream_analysis", "DEA"))
    params:
        gene_counts_dir = os.path.join(config["output_dir"], "gene_counts")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "DEA", "DEA.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "DEA", "DEA.txt")
    conda:
        "../envs/downstream.yaml"
    shell:
        """
        Rscript {input.r_script} \
            --conditionMapping {input.condition_mapping} \
            --geneCountsDir {params.gene_counts_dir} \
            --outputDir {output.output_dir}
        """

# rule gene_clustering:
#     """
    
#     """
#     input:
#         gene_counts = expand(rules.htseq_count.output.gene_counts, sample = SAMPLES),
#         condition_mapping = "conditionMapping.tsv",
#         r_script = os.path.join(config["scripts_dir"], "gene_clustering.R")
#     output:
#         flag_file = touch(os.path.join(config["output_dir"], "downstream_analysis", "gene_clustering", "clustering_done"))
#     params:
#         gene_counts_dir = os.path.join(config["output_dir"], "gene_counts"),
#         output_dir = os.path.join(config["output_dir"], "downstream_analysis", "data_exploration")
#     threads: 1
#     log:
#         os.path.join(config["logs_dir"], "gene_clustering", "gene_clustering.log")
#     benchmark:
#         os.path.join(config["benchmarks_dir"], "gene_clustering", "gene_clustering.txt")
#     # container:
#     shell:
#         """
#         Rscript {input.r_script} \
#             --conditionMapping {input.condition_mapping} \
#             --geneCountsDir {params.gene_counts_dir} \
#             --outputDir {params.output_dir}
#         """