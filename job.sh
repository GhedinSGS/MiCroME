#!/bin/bash

module load snakemake
snakemake -R --until FMAP_create_counts_matrix -j200 --rerun-incomplete --use-envmodules --latency-wait 300 --cluster "qsub -pe threaded {threads} -l h_vmem=10G -wd /hpcdata/lpd_sg/mchung/pipelines/SGSlab_metagenomics/stderrout/"

# snakemake -R --until spades -j200 --rerun-incomplete --use-envmodules --latency-wait 300 --cluster "qsub -pe threaded {threads} -l h_vmem=32G,himem -wd /hpcdata/lpd_sg/mchung/pipelines/SGSlab_metagenomics/stderrout/"

