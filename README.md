# Snakemake pipeline for the tumor macrophage project

This repository contains the Snakmake pipeline for the tumor macrophage project published in: DOI.

## Dependencies

Before running this pipeline, please make sure you have the following dependencies installed:
- Conda
- Singularity

The following section shows how to set up Snakemake using conda. In addition, all other dependencies will be handled by Snakemake using either conda or singularity. 

## Setting up Snakemake

You can set up the same version of Snakemake used in the analysis by using the `YAML` file in `envs/snakemake.yaml`. With conda installed, you can run:

```
conda env create -n snakemake -f envs/snakemake.yaml
```

## Configuring the pipeline

You can easily configure the pipeline using the `config.yaml` file at the root directory of this repository. The most important variable there is `raw_reads_dir`, which should be the path to a directory containing the forward and reverse FASTQ files. 

## Running the pipeline

Once set up, you can run the pipeline with:

```
snakemake -p -c <CORES> --rerun-triggers mtime --sdm conda apptainer
```

where `CORES` is the number of cores you want to run jobs in parallel.