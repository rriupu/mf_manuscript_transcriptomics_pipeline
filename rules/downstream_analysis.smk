rule data_exploration:
    """
    Run an initial data exploratory analysis with R.
    """
    input:
        gene_counts = expand(rules.htseq_count.output.gene_counts, sample = SAMPLES),
        condition_mapping = "conditionMapping.tsv",
        r_script = os.path.join("scripts", "data_exploration.R")
    output:
        flag_file = touch(os.path.join("results", "downstream_analysis", "data_exploration", "exploration_done"))
    params:
        gene_counts_dir = os.path.join("results", "gene_counts"),
        output_dir = os.path.join("results", "downstream_analysis", "data_exploration")
    threads: 1
    log:
        "logs/data_exploration/data_exploration.log"
    benchmark:
        "benchmarks/data_exploration/data_exploration.txt"
    # container:
    shell:
        """
        Rscript {input.r_script} \
            --conditionMapping {input.condition_mapping} \
            --geneCountsDir {params.gene_counts_dir} \
            --outputDir {params.output_dir}
        """

rule DEA:
    """
    Run a differential expression analysis with R.
    """
    input:
        gene_counts = expand(rules.htseq_count.output.gene_counts, sample = SAMPLES),
        condition_mapping = "conditionMapping.tsv",
        r_script = os.path.join("scripts", "DEA.R")
    output:
        flag_file = touch(os.path.join("results", "downstream_analysis", "DEA", "DEA_done"))
    params:
        gene_counts_dir = os.path.join("results", "gene_counts"),
        output_dir = os.path.join("results", "downstream_analysis", "DEA")
    threads: 1
    log:
        "logs/DEA/DEA.log"
    benchmark:
        "benchmarks/DEA/DEA.txt"
    # container:
    shell:
        """
        Rscript {input.r_script} \
            --conditionMapping {input.condition_mapping} \
            --geneCountsDir {params.gene_counts_dir} \
            --outputDir {params.output_dir}
        """

rule gene_clustering:
    """
    
    """
    input:
        gene_counts = expand(rules.htseq_count.output.gene_counts, sample = SAMPLES),
        condition_mapping = "conditionMapping.tsv",
        r_script = os.path.join("scripts", "gene_clustering.R")
    output:
        flag_file = touch(os.path.join("results", "downstream_analysis", "gene_clustering", "clustering_done"))
    params:
        gene_counts_dir = os.path.join("results", "gene_counts"),
        output_dir = os.path.join("results", "downstream_analysis", "data_exploration")
    threads: 1
    log:
        "logs/gene_clustering/gene_clustering.log"
    benchmark:
        "benchmarks/gene_clustering/gene_clustering.txt"
    # container:
    shell:
        """
        Rscript {input.r_script} \
            --conditionMapping {input.condition_mapping} \
            --geneCountsDir {params.gene_counts_dir} \
            --outputDir {params.output_dir}
        """

# rule TFBS_enrichment:
#     """
#     Run a differential expression analysis with R.
#     """
#     input:
#         gene_counts = expand(rules.htseq_count.output.gene_counts, sample = SAMPLES),
#         condition_mapping = "conditionMapping.tsv",
#         r_script = os.path.join("scripts", "data_exploration.R")
#     output:
#         touch(os.path.join("results", "downstream_analysis", "data_exploration", "exploration_done"))
#     params:
#         gene_counts_dir = os.path.join("results", "gene_counts"),
#         output_dir = os.path.join("results", "downstream_analysis", "data_exploration")
#     threads: 1
#     log:
#         "logs/{sample}/{sample}_DEA.log"
#     benchmark:
#         "benchmarks/{sample}/{sample}_DEA.txt"
#     container:
#     shell:
#         """
        
#         """