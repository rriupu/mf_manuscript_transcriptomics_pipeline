rule get_unibind_TFBS_enrichment_repo:
    """
    Clone the git repository with the UniBind TFBS enrichment tool.
    """
    output:
        unibind_enrichment_script = os.path.join(config["scripts_dir"], "unibind_enrichment", "bin", "UniBind_enrich.sh")
    params:
        clone_dir = os.path.join(config["scripts_dir"])
    threads: 1
    log:
        os.path.join(config["logs_dir"], "unibind_TFBS_enrichment", "get_repo.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "unibind_TFBS_enrichment", "get_repo.txt")
    shell:
        """
        cd {params.clone_dir}
        git clone https://bitbucket.org/CBGR/unibind_enrichment.git
        """

rule get_LOLA_database_for_unibind_enrichment:
    """
    """
    output:
        lola_db = os.path.join(config["data_dir"], "unibind_TFBS_enrichment", "hg38_robust_UniBind_LOLA.RDS")
    params:
        zenodo_fetch_script = os.path.join(config["scripts_dir"], "unibind_enrichment", "bin", "zenodo_fetch")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "unibind_TFBS_enrichment", "get_LOLA_db.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "unibind_TFBS_enrichment", "get_LOLA_db.txt")
    shell:
        """
        {params.zenodo_fetch_script} \
            4704641 \
            -f hg38_robust_UniBind_LOLA.RDS
        """

rule get_chrom_sizes:
    """
    """
    input:
        genome_fasta = rules.genome_fasta.output.genome_fasta
    output:
        chrom_sizes = os.path.join(config["data_dir"], "genome", "hg38.chromsizes")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "unibind_TFBS_enrichment", "get_chrom_sizes.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "unibind_TFBS_enrichment", "get_chrom_sizes.txt")
    container:
        "docker://staphb/samtools:1.22.1"
    shell:
        """
        samtools faidx {input.genome_fasta}
        cut -f1,2 {input.genome_fasta}.fai > {output.chrom_sizes}
        """

rule prepare_TSS_annotations:
    """
    """
    input:
        raw_TSS_annotations = "raw_TSS_annotations.bed",
        chrom_sizes = rules.get_chrom_sizes.output.chrom_sizes
    output:
        TSS_annotations = os.path.join(config["data_dir"], "unibind_TFBS_enrichment", "TSS_annotations", "mRNA_entrez_annotations_processed_sorted.bed")
    params:
        prepare_TSS_annotations_script = os.path.join(config["scripts_dir"], "prepare_TSS_annotations.sh"),
        output_dir = os.path.join(config["data_dir"], "unibind_TFBS_enrichment", "TSS_annotations")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "unibind_TFBS_enrichment", "prepare_TSS_annotations.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "unibind_TFBS_enrichment", "prepare_TSS_annotations.txt")
    conda:
        "../envs/downstream.yaml"
    shell:
        """
        bash {params.prepare_TSS_annotations_script} \
            {input.raw_TSS_annotations} \
            {input.chrom_sizes} \
            {params.output_dir}
        """

rule prepare_symbol2entrez_mapping:
    """
    """
    input:
        DEGs = os.path.join(config["output_dir"], "downstream_analysis", "DEA", "ADO_vs_M1", "DEA_results.tsv"),
    output:
        mapping = os.path.join(config["data_dir"], "unibind_TFBS_enrichment", "symbol2entrez_mapping.tsv")
    params:
        mapping_script = os.path.join(config["scripts_dir"], "generate_symbol2entrez_mapping.R")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "unibind_TFBS_enrichment", "prepare_symbol_to_entrez_mapping.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "unibind_TFBS_enrichment", "prepare_symbol_to_entrez_mapping.txt")
    shell:
        """
        Rscript {params.mapping_script} \
            -i {input.DEGs} \
            -o {output.mapping}
        """

rule generate_background:
    """
    """
    input:
        DEGs = os.path.join(config["output_dir"], "downstream_analysis", "DEA", "ADO_vs_M1", "DEA_results.tsv"),
        TSS_annotations = rules.prepare_TSS_annotations.output.TSS_annotations,
        symbol2entrez_mapping = rules.prepare_symbol2entrez_mapping.output.mapping
    output:
        background_file = os.path.join(config["output_dir"], "unibind_enrichment", "background.bed")
    params:
        output_dir = os.path.join(config["output_dir"], "unibind_enrichment"),
        generate_background_script = os.path.join(config["scripts_dir"], "prepare_background.R")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "unibind_TFBS_enrichment", "generate_background.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "unibind_TFBS_enrichment", "generate_background.txt")
    shell:
        """
        Rscript {params.generate_background_script} \
            -i {input.DEGs} \
            -b {input.TSS_annotations} \
            -m {input.symbol2entrez_mapping} \
            -o {params.output_dir}
        """

rule generate_foreground:
    """
    """
    input:
        DEGs = os.path.join(config["output_dir"], "downstream_analysis", "DEA", "{comparison}", "DEA_results.tsv"),
        background_file = rules.generate_background.output.background_file,
        symbol2entrez_mapping = rules.prepare_symbol2entrez_mapping.output.mapping
    output:
        expand(os.path.join(config["output_dir"], "unibind_enrichment", "{comparison}", "foreground_{deg_subset}.bed"), deg_subset = ["all", "upregulated", "downregulated"])
    params:
        generate_foreground_script = os.path.join(config["scripts_dir"], "prepare_foreground.R"),
        output_dir = os.path.join(config["output_dir"], "unibind_enrichment", "{comparison}")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "unibind_TFBS_enrichment", "generate_foreground_{comparison}.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "unibind_TFBS_enrichment", "generate_foreground_{comparison}.txt")
    shell:
        """
        Rscript {params.generate_foreground_script} \
            -i {input.DEGs} \
            -b {input.background_file} \
            -m {input.symbol2entrez_mapping} \
            -o {params.output_dir}
        """

rule TFBS_enrichment:
    """
    Run the UniBind TFBS enrichment analysis for each comparison's DEGs.
    """
    input:
        enrichment_script = rules.get_unibind_TFBS_enrichment_repo.output.unibind_enrichment_script,
        lola_db = rules.get_LOLA_database_for_unibind_enrichment.output.lola_db,
        foreground_file = rules.generate_foreground.output.foreground_file,
        background_file = rules.generate_background.output.background_file,
        DEA = get_comparisons
    output:
        enrichment_results = os.path.join(config["output_dir"], "unibind_enrichment", "{comparison}", "{deg_subset}", "allEnrichments.tsv")
    params:
        output_dir = os.path.join(config["output_dir"], "unibind_enrichment", "{comparison}", "{deg_subset}")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "unibind_TFBS_enrichment", "{comparison}_{deg_subset}_enrichment.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "unibind_TFBS_enrichment", "{comparison}_{deg_subset}_enrichment.txt")
    container:
        "docker://cbgr/unibind_enrichment:0.2"
    shell:
        """
        {params.enrichment_script} \
            oneSetBg \
            {input.lola_db} \
            {input.foreground_file} \
            {input.background_file} \
            {params.output_dir}
        """

def aggregate_enrichment_results(wildcards):
    checkpoint_output = checkpoints.DEA.get(**wildcards).output[0]
    out = expand(
        os.path.join(
            config["output_dir"],
            "unibind_enrichment",
            "{comparison}",
            "{deg_subset}",
            "allEnrichments.tsv"
        ), comparison = glob_wildcards(
            os.path.join(
                checkpoint_output,
                "{comparison}",
                "DEA_results.tsv")).comparison,
        , deg_subset = ["all", "upregulated", "downregulated"]
    )
    return out

rule generate_TFBS_enrichment_summary_plots:
    """
    """
    input:
        aggregate_enrichment_results
    output:
        summary_plots = expand(os.path.join(config["output_dir"], "unibind_enrichment", "summary_plots", "{deg_subset}", "summary_dotplot_all.pdf"), deg_subset = ["all", "upregulated", "downregulated"]),
        upset_plots = expand(os.path.join(config["output_dir"], "unibind_enrichment", "upset_plots", "{deg_subset}", "upset_plot_all_comparisons.pdf"), deg_subset = ["all", "upregulated", "downregulated"])
    params:
        TFBS_enrichment_summary_plot = os.path.join(config["scripts_dir"], "TFBS_enrichment_summary.R"),
        input_dir = os.path.join(config["output_dir"], "unibind_enrichment", "{deg_subset}"),
        summary_output_dir = os.path.join(config["output_dir"], "unibind_enrichment", "summary_plots"),
        upset_plots_output_dir = os.path.join(config["output_dir"], "unibind_enrichment", "summary_plots")
    threads: 1
    log:
        os.path.join(config["logs_dir"], "unibind_TFBS_enrichment", "TFBS_summary_plots_{comparison}_{deg_subset}.log")
    benchmark:
        os.path.join(config["benchmarks_dir"], "unibind_TFBS_enrichment", "TFBS_summary_plots_{comparison}_{deg_subset}.txt")
    shell:
        """
        for comparison in all downregulated upregulated; do
            
            Rscript {params.TFBS_enrichment_summary_plot} \
                -i {params.input_dir} \
                -c $comparison \
                -o {params.summary_output_dir}/$comparison/

            Rscript {params.TFBS_enrichment_upset_plot} \
                -i {params.input_dir} \
                -c $comparison \
                -o {params.upset_plots_output_dir}/$comparison/

        done
        """