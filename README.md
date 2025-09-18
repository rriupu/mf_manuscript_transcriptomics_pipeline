# Snakemake pipeline for the tumor macrophage project

This repository contains the Snakmake pipeline for the tumor macrophage project published in: DOI.

## Dependencies

Conda
Singularity

## Setting up Snakemake

You can set up the same version of Snakemake used in the analysis by using the `YAML` file in `envs/snakemake.yaml`. With conda installed, you can run:

```
conda env create -n snakemake -f envs/snakemake.yaml
```