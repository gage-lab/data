# Snakemake workflow: ngs-test-data

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.16-brightgreen.svg)](https://snakemake.readthedocs.org/)
[![Tests](https://github.com/gage-lab/ngs-test-data/actions/workflows/main.yaml/badge.svg)](https://github.com/gage-lab/ngs-test-data/actions/workflows/main.yaml)

This workflow creates small test datasets for NGS data analyses.

## Authors

- Michael Cuoco (@mikecuoco), https://michaelcuoco.com

## Usage

Use this repo as a submodule in your project! To generate the test data, simply run the following from this repo's base

```bash
# generate bulk RNA-seq test data
snakemake rnaseq --use-conda -c1

# generate 10x v3 scRNA-seq test data
snakemake scrnaseq_10x_v3 --use-conda -c1

# generate Whole-Genome Sequencing test data
snakemake wgs --use-conda -c1
```
