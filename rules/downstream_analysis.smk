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
        r_script = os.path.join("scripts", "data_exploration.R"),
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
            --outputDir {params.output_dir} \
        2> {log}
        """

rule gene_clustering:
    """
    
    """
    input:
        gene_counts = expand(rules.htseq_count.output.gene_counts, sample = SAMPLES),
        condition_mapping = "conditionMapping.tsv"
    output:
        clusters = os.path.join(config["output_dir"], "downstream_analysis", "gene_clustering", "clusters.pdf")
    params:
        r_script = os.path.join("scripts", "gene_clustering.R"),
        gene_counts_dir = os.path.join(config["output_dir"], "gene_counts"),
        output_dir = os.path.join(config["output_dir"], "downstream_analysis", "gene_clustering")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "gene_clustering", "gene_clustering.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "gene_clustering", "gene_clustering.txt")
    conda:
        "../envs/downstream.yaml"
    shell:
        """
        Rscript {params.r_script} \
            --conditionMapping {input.condition_mapping} \
            --geneCountsDir {params.gene_counts_dir} \
            --outputDir {params.output_dir} \
        2> {log}
        """

rule GSVA:
    """
    """
    input:
        gene_counts = expand(rules.htseq_count.output.gene_counts, sample = SAMPLES),
        condition_mapping = "conditionMapping.tsv"
    output:
        GSVA_plot = os.path.join(config["output_dir"], "downstream_analysis", "GSVA", "GSVA_GO.pdf")
    params:
        gsva_script = os.path.join("scripts", "GSVA.R"),
        gene_counts_dir = os.path.join(config["output_dir"], "gene_counts"),
        output_dir = os.path.join(config["output_dir"], "downstream_analysis", "GSVA")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "GSVA", "GSVA.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "GSVA", "GSVA.txt")
    conda:
        "../envs/downstream.yaml"
    shell:
        """
        Rscript {params.gsva_script} \
            -c {input.condition_mapping} \
            -g {params.gene_counts_dir} \
            -o {params.output_dir} \
        2> {log}
        """

checkpoint DEA:
    """
    Run a differential expression analysis with R.
    """
    input:
        gene_counts = expand(rules.htseq_count.output.gene_counts, sample = SAMPLES),
        condition_mapping = "conditionMapping.tsv",
        r_script = os.path.join("scripts", "DEA.R")
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
            --outputDir {output.comparison_results} \
        2> {log}
        """