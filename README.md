# Snakemake workflow: ngs-test-data

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.16-brightgreen.svg)](https://snakemake.readthedocs.org/)

This workflow creates small test datasets for NGS data analyses.

## Authors

* Johannes Köster (@johanneskoester), https://koesterlab.github.io
* Michael Cuoco (@mikecuoco), https://michaelcuoco.com

## Usage

Use this repo as a submodule in your project! To generate the test data, simply run the following from this repo's base

```bash
snakemake all --use-conda -c1
```
