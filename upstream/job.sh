#!/bin/bash

snakemake --unlock
snakemake -j200 --rerun-incomplete --use-envmodules --use-conda --latency-wait 300 --cluster "sbatch -J sasha_mushegian_12 --time=24:00:00 --mem-per-cpu=36000 --cpus-per-task={threads} --chdir /data/lpd_sg/mchung/pipelines/SGSlab_metagenomics/stderrout/"