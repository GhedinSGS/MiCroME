#!/bin/bash

snakemake --unlock
snakemake  -j200 --rerun-incomplete --use-envmodules --use-conda --latency-wait 300 --cluster "sbatch -J christina_yek_02 --time=24:00:00 --mem-per-cpu=36000 --cpus-per-task={threads} --chdir /data/lpd_sg/mchung/pipelines/SGSlab_metagenomics/stderrout/"

# snakemake -R --until spades -j200 --rerun-incomplete --use-envmodules --latency-wait 300 --cluster "qsub -pe threaded {threads} -l h_vmem=32G,himem -wd /hpcdata/lpd_sg/mchung/pipelines/SGSlab_metagenomics/stderrout/"

