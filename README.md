# Snakemake workflow: ngs-test-data

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.16-brightgreen.svg)](https://snakemake.readthedocs.org/)
[![Tests](https://github.com/gage-lab/ngs-test-data/actions/workflows/main.yml/badge.svg)](https://github.com/gage-lab/ngs-test-data/actions/workflows/main.yml)

This workflow creates small test datasets for NGS data analyses.

## Authors

- Johannes Köster (@johanneskoester), https://koesterlab.github.io
- Michael Cuoco (@mikecuoco), https://michaelcuoco.com

## Usage

Use this repo as a submodule in your project! To generate the test data, simply run the following from this repo's base

```bash
# generate rnaseq test data
snakemake rnaseq --use-conda -c1
```
