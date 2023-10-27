# SGS Lab Metagenomics Analysis Pipeline
# Usage: snakemake -j200 -c4 --rerun-incomplete --use-envmodules --latency-wait 300 --cluster "qsub -wd /hpcdata/lpd_sg/mchung/pipelines/SGSlab_metagenomics/stderrout/"
#
# Author: Matthew Chung < chungm6@nih.gov >
# Laboratory of Parasitic Diseases, Systems Genomics Section, NIAID, NIH

configfile: "config.yaml"
version: "1.00"

import pandas as pd

# SAMPLE, = glob_wildcards(config["params"]["reads_dir"] + "/{ids}.1.fastq.gz")
SAMPLE, = glob_wildcards(config["params"]["reads_dir"] + "/{ids}_R1_001.fastq.gz")
# SAMPLE, = glob_wildcards(config["params"]["reads_dir"] + "/{ids}_R1_001.trimmed.fastq.gz")

SAMPLE = [x for x in SAMPLE if "Undetermined" not in x]

# print(SAMPLE)

wildcard_constraints:
    sample='|'.join([re.escape(x) for x in SAMPLE]),

rule all:
    input:
        config["params"]["output"] + "tables/counts.ko.ortho.fmap.tsv",
        config["params"]["output"] + "tables/counts.ko.module.fmap.tsv",
        config["params"]["output"] + "tables/counts.ko.pathway.fmap.tsv",
        config["params"]["output"] + "tables/counts.taxa.bracken.tsv",
        config["params"]["output"] + "tables/annotation.taxa.bracken.tsv",
        config["params"]["output"] + "tables/counts.coverM.MAG.tsv",
        config["params"]["output"] + "tables/counts.coverM.viralMAG.tsv",
        config["params"]["output"] + "tables/counts.coverM.CARD.tsv",
        config["params"]["output"] + "tables/counts.coverM.MEGARes.tsv",
        config["params"]["output"] + "tables/MAG_summary.tsv",
        config["params"]["output"] + "ResMiCo/report.txt",
        expand(config["params"]["output"] + "MetaPhlAn/{sample}.MetaPhlAn.txt",sample=SAMPLE),
        expand(config["params"]["output"] + "MetaCRAST/{sample}/totalSpacers.fa",sample=SAMPLE),
        expand(config["params"]["output"] + "HUMAnN/{sample}/test.txt",sample=SAMPLE)

rule fastp:
    input:
        fq1 = config["params"]["reads_dir"] + "{sample}_R1_001.fastq.gz",
        fq2 = config["params"]["reads_dir"] + "{sample}_R2_001.fastq.gz"
    output:
        fq1 = config["params"]["output"] + "fastp/{sample}.1.fq.gz",
        fq2 = config["params"]["output"] + "fastp/{sample}.2.fq.gz",
        qcfile = config["params"]["output"] + "fastp/{sample}.html"
    threads: 8
    envmodules:
        "fastp/0.23.4"
    shell:
        """
        fastp --dedup \
            -i {input.fq1} -I {input.fq2} \
            -o {output.fq1} -O {output.fq2} \
            -h {output.qcfile} -w {threads} \
        """

rule Bowtie2_deplete_host_reads:
    input:
        fq1 = config["params"]["output"] + "fastp/{sample}.1.fq.gz",
        fq2 = config["params"]["output"] + "fastp/{sample}.2.fq.gz"
    params:
        dir = config["params"]["output"] + "Bowtie2/",
        ref = config["ref"]["host_ref"]
    output:
        fq1 = config["params"]["output"] + "Bowtie2/{sample}.exclude.host.1.fq.gz",
        fq2 = config["params"]["output"] + "Bowtie2/{sample}.exclude.host.2.fq.gz"
    threads: 8
    envmodules:
        "bowtie2/2.4.2-Python-3.6.12"
    shell:
        """
        bowtie2 \
            -p {threads} \
            -x {params.ref} \
            -1 {input.fq1} \
            -2 {input.fq2} \
            --un-conc-gz {params.dir}/{wildcards.sample} \
            -S {params.dir}/{wildcards.sample}.sam

        rm {params.dir}/{wildcards.sample}.sam

        mv {params.dir}/{wildcards.sample}.1 {output.fq1}
        mv {params.dir}/{wildcards.sample}.2 {output.fq2}
        """

rule Kraken2:
    input:
        fq1 = config["params"]["output"] + "Bowtie2/{sample}.exclude.host.1.fq.gz",
        fq2 = config["params"]["output"] + "Bowtie2/{sample}.exclude.host.2.fq.gz"
    params:
        db = config["ref"]["kraken2_db"]
    output:
        kraken2_report = config["params"]["output"] + "Kraken2/{sample}.report.txt",
        kraken2_output = config["params"]["output"] + "Kraken2/{sample}.output.tsv"
    threads: 8
    envmodules:
        "kraken/2.1.2-gompi-2021a"
    shell:
        """
        if [ $(stat -c %s {input.fq1}) -gt 20 ]
        then
            kraken2 --threads {threads} --db {params.db} \
                --report {output.kraken2_report} --output {output.kraken2_output} \
                --paired --gzip-compressed {input.fq1} {input.fq2}
        else
            touch {output.kraken2_report}
            touch {output.kraken2_output}
        fi
        """

rule Bracken:
    input:
        kraken2_report = config["params"]["output"] + "Kraken2/{sample}.report.txt"
    params:
        db = config["ref"]["kraken2_db"]
    output:
        bracken_output = config["params"]["output"] + "Bracken/{sample}.output.txt"
    envmodules:
        "bracken/2.5-foss-2016b"
    shell:
        """
        if [ $(stat -c %s {input.kraken2_report}) -gt 0 ]
        then
            bracken -d {params.db} -i {input.kraken2_report} -o {output.bracken_output} -r 100 -l S
        else
            touch {output.bracken_output}
        fi
        """

rule Bracken_create_counts_matrix:
    input:
        tsv = expand(config["params"]["output"] + "Bracken/{sample}.output.txt",sample=SAMPLE)
    params:
        bin = config["bin"]["scripts"],
        dir = config["params"]["output"]
    output:
        tsv1 = config["params"]["output"] + "tables/counts.taxa.bracken.tsv",
        tsv2 = config["params"]["output"] + "tables/annotation.taxa.bracken.tsv"
    threads: 1
    envmodules:
        "R/4.2.1"
    shell:
        """
        Rscript {params.bin}/bracken_create_counts_matrix.R -d {params.dir}
        """

rule deplete_human_reads:
    input:
        fq1 = config["params"]["output"] + "Bowtie2/{sample}.exclude.host.1.fq.gz",
        fq2 = config["params"]["output"] + "Bowtie2/{sample}.exclude.host.2.fq.gz",
        kraken2_output = config["params"]["output"] + "Kraken2/{sample}.output.tsv",
        kraken2_report = config["params"]["output"] + "Kraken2/{sample}.report.txt"
    params:
        krakentools = config["bin"]["krakentools"],
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq"
    output:
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq.gz",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq.gz"
    envmodules:
        "python/3.7.3-foss-2016b"
    shell:
        """
        if [ $(stat -c %s {input.fq1}) -gt 20 ]
        then
            rm -f {params.fq1}
            rm -f {params.fq2}

            zcat {input.fq1} | split -l 120000000 -d - {input.fq1}
            zcat {input.fq2} | split -l 120000000 -d - {input.fq2}

            cat {input.kraken2_output} | split -l 20000000 -d - {input.kraken2_output}

            for x in $(ls {input.fq1}[0-9]* | sed "s/.*.1.fq.gz//g" | sort -n)
            do
                for y in $(ls {input.kraken2_output}[0-9]* | sed "s/.*output.tsv//g"| sort -n)
                do
                    {params.krakentools}/extract_kraken_reads.py --exclude --include-children --fastq-output \
                        -t 9606 \
                        --max 1000000000 \
                        -s1 {input.fq1}"$x" \
                        -s2 {input.fq2}"$x" \
                        -k {input.kraken2_output}"$y" \
                        -r {input.kraken2_report} \
                        -o {params.fq1}"$x"_"$y" \
                        -o2 {params.fq2}"$x"_"$y"

                    cat {params.fq1}"$x"_"$y" >> {params.fq1}
                    cat {params.fq2}"$x"_"$y" >> {params.fq2}

                    rm {params.fq1}"$x"_"$y"
                    rm {params.fq2}"$x"_"$y"
                done
                rm {input.fq1}"$x"
                rm {input.fq2}"$x"
            done

            rm {input.kraken2_output}[0-9]*

            gzip {params.fq1}
            gzip {params.fq2}
        else
            touch {output.fq1}
            touch {output.fq2}
        fi
        """

rule MEGARes_coverM:
    input:
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq.gz",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq.gz"
    params:
        db = config["ref"]["megares_db"]
    output:
        tsv = config["params"]["output"] + "coverM/MEGARes/{sample}.tsv"
    threads: 8
    envmodules:
        "coverm/0.6.1",
        "samtools/1.9-goolf-1.7.20"
    shell:
        """
        coverm contig \
            --methods count \
            --coupled {input.fq1} {input.fq2} \
            --reference {params.db} \
            -t {threads} \
            -o {output.tsv} \
            --exclude-supplementary  
        """

rule CARD_coverM:
    input:
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq.gz",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq.gz"
    params:
        db = config["ref"]["card_db"]
    output:
        tsv = config["params"]["output"] + "coverM/CARD/{sample}.tsv"
    threads: 8
    envmodules:
        "coverm/0.6.1",
        "samtools/1.9-goolf-1.7.20"
    shell:
        """
        coverm contig \
            --methods count \
            --coupled {input.fq1} {input.fq2} \
            --reference {params.db} \
            -t {threads} \
            -o {output.tsv} \
            --exclude-supplementary  
        """

rule MEGARes_coverM_create_counts_matrix:
    input:
        tsv = expand(config["params"]["output"] + "coverM/MEGARes/{sample}.tsv",sample=SAMPLE)
    params:
        bin = config["bin"]["scripts"],
        dir = config["params"]["output"]
    output:
        tsv = config["params"]["output"] + "tables/counts.coverM.MEGARes.tsv"
    envmodules:
        "R/4.2.1"
    shell:
        """
        Rscript {params.bin}/coverM_MEGARes_create_counts_matrix.R -d {params.dir}
        """

rule CARD_coverM_create_counts_matrix:
    input:
        tsv = expand(config["params"]["output"] + "coverM/CARD/{sample}.tsv",sample=SAMPLE)
    params:
        bin = config["bin"]["scripts"],
        dir = config["params"]["output"]
    output:
        tsv = config["params"]["output"] + "tables/counts.coverM.CARD.tsv"
    envmodules:
        "R/4.2.1"
    shell:
        """
        Rscript {params.bin}/coverM_CARD_create_counts_matrix.R -d {params.dir}
        """

rule AMRFinderPlus:
    input:
        fa = config["params"]["output"] + "SPAdes/{sample}/contigs.final.fasta"
    output:
        tsv = config["params"]["output"] + "AMRFinderPlus/{sample}.tsv"
    threads: 8
    envmodules:
        "amrfinderplus/3.11.11"
    shell:
        """
        amrfinder \
            --nucleotide {input.fa} \
            --output {output.tsv} \
            --threads {threads}
        """

rule FMAP_mapping:
    input:
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq.gz",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq.gz"
    params:
        fmap = config["bin"]["fmap"]
    output:
        blastx = config["params"]["output"] + "FMAP/{sample}.blastx_hits.txt"
    threads: 8
    envmodules:
        "fmap/0.15-goolf-1.7.20"
    shell:
        """
        perl {params.fmap}/FMAP_mapping.pl \
            -p {threads} \
            {input.fq1} {input.fq2} > {output.blastx}
        """

rule FMAP_quantification:
    input:
        blastx = config["params"]["output"] + "FMAP/{sample}.blastx_hits.txt"
    params:
        fmap = config["bin"]["fmap"]
    output:
        abundance = config["params"]["output"] + "FMAP/{sample}.abundance.txt"
    threads: 1
    envmodules:
        "fmap/0.15-goolf-1.7.20"
    shell:
        """
        perl {params.fmap}/FMAP_quantification.pl {input.blastx} > {output.abundance}
        """

rule FMAP_create_counts_matrix:
    input:
        tsv = expand(config["params"]["output"] + "FMAP/{sample}.abundance.txt",sample=SAMPLE)
    params:
        bin = config["bin"]["scripts"],
        dir1 = config["params"]["output"],
        dir2 = config["bin"]["fmap"]
    output:
        tsv1 = config["params"]["output"] + "tables/counts.ko.ortho.fmap.tsv",
        tsv2 = config["params"]["output"] + "tables/counts.ko.module.fmap.tsv",
        tsv3 = config["params"]["output"] + "tables/counts.ko.pathway.fmap.tsv"
    threads: 1
    envmodules:
        "R/4.2.1"
    shell:
        """
        Rscript {params.bin}/fmap_create_counts_matrix.R -d {params.dir1} -f {params.dir2}
        """

rule SPAdes:
    input:
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq.gz",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq.gz"
    params:
        dir = config["params"]["output"] + "SPAdes/{sample}/"
    output:
        fa = config["params"]["output"] + "SPAdes/{sample}/contigs.fasta"
    threads: 8
    envmodules:
        "spades/3.15.2"
    shell:
    	"""
        if [ -s {input.fq1} ]
        then
            spades.py --meta --only-assembler -t {threads} -1 {input.fq1} -2 {input.fq2} -o {params.dir}

            sed -i "s/^>/>{wildcards.sample}|/g" {output.fa}
        else
            mkdir -p {params.dir}
            touch {output.fa}
        fi
        """

rule Kraken2_MAG:
    input:
        fa = config["params"]["output"] + "SPAdes/{sample}/contigs.fasta"
    params:
        db = config["ref"]["kraken2_db"],
    output:
        kraken2_report = config["params"]["output"] + "Kraken2_MAG/{sample}.report.txt",
        kraken2_output = config["params"]["output"] + "Kraken2_MAG/{sample}.output.tsv"
    threads: 8
    envmodules:
        "kraken/2.1.2-gompi-2021a"
    shell:
        """
        if [ -s {input.fa} ]
        then
            kraken2 --threads {threads} --db {params.db} \
                --report {output.kraken2_report} --output {output.kraken2_output} \
                {input.fa}
        else
            touch {output.kraken2_report}
            touch {output.kraken2_output}
        fi
        """

rule deplete_human_MAGs:
    input:
        fa = config["params"]["output"] + "SPAdes/{sample}/contigs.fasta",
        kraken2_output = config["params"]["output"] + "Kraken2_MAG/{sample}.output.tsv"
    params:
        bin = config["bin"]["scripts"],
        bp = config["params"]["min_MAG_bp_size"],
        dir = config["params"]["output"] + "SPAdes/{sample}/"
    output:
        fa = config["params"]["output"] + "SPAdes/{sample}/contigs.final.fasta",
    envmodules:
        "samtools/1.9-goolf-1.7.20"
    shell:
        """
        if [ -s {input.fa} ]
        then
            cat {input.kraken2_output} | awk '$3 != 9606 {{print $2}}' > {params.dir}/nonhuman_MAGs.tsv
            if [ -s {params.dir}/nonhuman_MAGs.tsv ]
            then
                samtools faidx {input.fa} -r {params.dir}/nonhuman_MAGs.tsv > {params.dir}/contigs.nonhuman.fasta
            else
                touch {params.dir}/contigs.nonhuman.fasta
            fi

            if [ -s {params.dir}/contigs.nonhuman.fasta ]
            then
                perl {params.bin}/residues.pl {params.dir}/contigs.nonhuman.fasta | awk '$2 > {params.bp} {{print $1}}' > {params.dir}/gt{params.bp}bp_MAGs.tsv
            else
                touch {params.dir}/gt{params.bp}bp_MAGs.tsv
            fi

            if [ -s {params.dir}/gt{params.bp}bp_MAGs.tsv ]
            then
                samtools faidx {input.fa} -r {params.dir}/gt{params.bp}bp_MAGs.tsv > {output.fa} 
            else
                touch {output.fa}
            fi
        else
            touch {output.fa}
        fi
        """

rule ResMiCo_map:
    input:
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq.gz",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq.gz",
        ref = config["params"]["output"] + "SPAdes/{sample}/contigs.final.fasta"
    output:
        bam = config["params"]["output"] + "ResMiCo/bam/{sample}.bam"
    threads: 4
    envmodules:
        "bowtie2/2.4.2-Python-3.6.12",
        "samtools/1.9-goolf-1.7.20"
    shell:
        """
        if [ -s {input.ref} ]
        then
            bowtie2-build --large-index --threads {threads} {input.ref} {input.ref}

            bowtie2 \
                -p {threads} \
                -x {input.ref} \
                -1 {input.fq1} \
                -2 {input.fq2} | samtools sort -o {output.bam} -
        else
            touch {output.bam}
        fi
        """

rule ResMiCo_evaluate:
    input:
        bam = expand(config["params"]["output"] + "ResMiCo/bam/{sample}.bam",sample=SAMPLE)
    params:
        dir1 = config["params"]["output"] + "SPAdes/",
        dir2 = config["params"]["output"] + "ResMiCo/"
    output:
        map = config["params"]["output"] + "ResMiCo/map.tsv",
        tsv = config["params"]["output"] + "ResMiCo/report.txt"
    threads: 8
    envmodules:
        "resmico/1.2.2"
    shell:
        """
        echo -e "Taxon\tFasta\tSample\tBAM" > {output.map}
        for BAM in $(find {params.dir2}/bam/ -name "*[.]bam")
        do
            SAMPLE=$(echo $BAM | sed "s/.*\\///g" | sed "s/[.]bam//g")
            echo -e "$SAMPLE\t{params.dir1}/$SAMPLE/contigs.final.fasta\t$SAMPLE\t$BAM" >> {output.map}
        done
        cat {output.map} | tr -s '[:blank:]' '\t' > {params.dir2}/temp
        mv {params.dir2}/temp {output.map}

        resmico bam2feat --n-threads {threads} --outdir {params.dir2}/features {output.map}

        resmico evaluate \
            --n-procs {threads} \
            --min-avg-coverage 1 \
            --save-path {params.dir2}/predictions \
            --save-name {params.dir2}/default_model \
            --feature-files-path {params.dir2}/features
        """

rule MetaPhlAn_MAG:
    input:
        fa = config["params"]["output"] + "SPAdes/{sample}/contigs.fasta"
    params:
        db = config["ref"]["metaphlan_db"]
    output:
        bowtie2out = config["params"]["output"] + "MetaPhlAn_MAG/{sample}.bowtie2out.txt",
        txt = config["params"]["output"] + "MetaPhlAn_MAG/{sample}.MetaPhlAn.txt"
    threads: 8
    envmodules:
        "metaphlan/4.0.2"
    shell:
        """
        metaphlan \
            {input.fa} \
            --bowtie2db {params.db} \
            --bowtie2out {output.bowtie2out} \
            -o {output.txt} \
            --nproc {threads} \
            -t rel_ab_w_read_stats \
            --add_viruses --input_type fasta
        """

rule VAMB_concantenate:
    input:
        fa = expand(config["params"]["output"] + "SPAdes/{sample}/contigs.final.fasta",sample=SAMPLE),
        kraken2_output = expand(config["params"]["output"] + "Kraken2_MAG/{sample}.output.tsv",sample=SAMPLE)
    params:
        bin = config["bin"]["scripts"],
        dir = config["params"]["output"] + "VAMB/"
    output:
        fa = config["params"]["output"] + "VAMB/catalogue.fna",
        tsv = config["params"]["output"] + "VAMB/catalogue.residues.tsv"
    envmodules:
        "samtools/1.9-goolf-1.7.20",
        "vamb/3.0.9"
    shell:
        """
        cat {input.fa} > {params.dir}/catalogue.fna
        perl {params.bin}/residues.pl {output.fa} > {output.tsv}
        """

rule VAMB_map:
    input:
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq.gz",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq.gz",
        fna = config["params"]["output"] + "VAMB/catalogue.fna"
    output:
        bam = config["params"]["output"] + "VAMB/bam/{sample}.bam"
    threads: 8
    envmodules:
        "minimap2/2.24",
        "samtools/1.9-goolf-1.7.20"
    shell:
        """
        minimap2 -t {threads} -N 5 -ax sr {input.fna} {input.fq1} {input.fq2} | \
            samtools view -h -F 3584 -@ {threads} -o {output.bam} -
        """

rule VAMB_cluster:
    input:
        bam = expand(config["params"]["output"] + "VAMB/bam/{sample}.bam",sample=SAMPLE),
        fna = config["params"]["output"] + "VAMB/catalogue.fna"
    params:
        bp = config["params"]["min_MAG_bp_size"],
        dir1 = config["params"]["output"] + "VAMB/bam",
        dir2 = config["params"]["output"] + "VAMB/clusters"
    output:
        tsv = config["params"]["output"] + "VAMB/clusters/clusters.tsv"
    envmodules:
        "vamb/3.0.9"
    shell:
        """
        rm -rf {params.dir2}
        vamb --outdir {params.dir2} --fasta {input.fna} --bamfiles {params.dir1}/*bam -o "___"  --minfasta {params.bp}
        """

rule GTDB_Tk:
    input:
        tsv = config["params"]["output"] + "VAMB/clusters/clusters.tsv"
    params:
        dir1 = config["params"]["output"] + "VAMB/clusters/bins/",
        dir2 = config["params"]["output"] + "GTDB-Tk/"
    output:
        tsv = config["params"]["output"] + "GTDB-Tk/classify/gtdbtk.bac120.summary.tsv"
    threads: 8
    envmodules:
        "gtdbtk/2.1.1"
    shell:
        """
        gtdbtk classify_wf --genome_dir {params.dir1} --out_dir {params.dir2} --cpus {threads}
        """

rule VirSorter2_pass1:
    input:
        fna=config["params"]["output"] + "SPAdes/{sample}/contigs.final.fasta"
    params:
        bin = config["bin"]["scripts"],
        dir=config["params"]["output"] + "VirSorter2/pass1/{sample}/",
        bp = config["params"]["min_MAG_bp_size"]
    output:
        fna=config["params"]["output"] + "VirSorter2/pass1/{sample}/final-viral-combined.fa"
    threads: 8
    envmodules:
        "virsorter/2.2.3-Python-3.8.10"
    shell:
        """
        contigs=$(perl {params.bin}/residues.pl {input.fna} | awk '$2 > {params.bp} {{print $1}}' | wc -l)

        if [[ $contigs -gt 0 ]]
        then
            virsorter run -i {input.fna} -w {params.dir} --include-groups dsDNAphage,ssDNA --keep-original-seq --min-length {params.bp} -j {threads} all
        else
            touch {output.fna}
        fi
        """

rule CheckV:
    input:
        fna=config["params"]["output"] + "VirSorter2/pass1/{sample}/final-viral-combined.fa"
    params:
        dir=config["params"]["output"] + "CheckV/{sample}/"
    output:
        fna=config["params"]["output"] + "CheckV/{sample}/combined_viruses.fna"
    threads: 8
    envmodules:
        "checkv/0.8.1"
    shell:
        """
        if [ -s {input.fna} ]
        then
            checkv end_to_end {input.fna} {params.dir} -t {threads}
            cat {params.dir}/proviruses.fna {params.dir}/viruses.fna > {output.fna}
        else
            touch {output.fna}
        fi
        """

rule VirSorter2_pass2:
    input:
        fna = config["params"]["output"] + "CheckV/{sample}/combined_viruses.fna"
    params:
        dir = config["params"]["output"] + "VirSorter2/pass2/{sample}/",
        bp = config["params"]["min_MAG_bp_size"]
    output:
        fna1 = config["params"]["output"] + "VirSorter2/pass2/{sample}/final-viral-combined.fa",
        fna2 = config["params"]["output"] + "VirSorter2/pass2/{sample}/for-dramv/final-viral-combined-for-dramv.fa",
        tab = config["params"]["output"] + "VirSorter2/pass2/{sample}/for-dramv/viral-affi-contigs-for-dramv.tab"
    threads: 8
    envmodules:
        "virsorter/2.2.3-Python-3.8.10"
    shell:
        """
        if [ -s {input.fna} ]
        then
            touch {output.fna2}
            virsorter run -i {input.fna} -w {params.dir} --include-groups dsDNAphage,ssDNA --seqname-suffix-off --viral-gene-enrich-off --prep-for-dramv --min-length {params.bp} -j {threads} all
        else
            mkdir -p {params.dir}/for-dramv/

            touch {output.fna1}
            touch {output.fna2}
            touch {output.tab}
        fi
        """

rule create_viralMAG_ref:
    input:
        fna = expand(config["params"]["output"] + "VirSorter2/pass2/{sample}/final-viral-combined.fa",sample=SAMPLE)
    params:
        dir = config["params"]["output"]
    output:
        fna = config["params"]["output"] + "VirSorter2/viral_MAGs.fna"
    shell:
        """
        cat {input.fna} > {output.fna}
        """

rule iPHoP:
    input:
        fna = config["params"]["output"] + "VirSorter2/pass2/{sample}/final-viral-combined.fa"
    params:
        dir = config["params"]["output"] + "iPHoP/{sample}/",
        db = config["ref"]["iphop_db"]
    output:
        csv = config["params"]["output"] + "iPHoP/{sample}/Host_prediction_to_genus_m90.csv"
    threads: 8
    envmodules:
        "iphop/1.3.0",
        "samtools/1.9-goolf-1.7.20"
    shell:
        """
        if [ -s {input.fna} ]
        then
            iphop predict --fa_file {input.fna} --out_dir {params.dir} --db_dir {params.db} -t {threads} || true
        fi
        touch {output.csv}
        """

rule coverM_MAG:
    input:
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq.gz",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq.gz",
        ref = config["params"]["output"] + "VAMB/catalogue.fna"
    output:
        tsv = config["params"]["output"] + "coverM/MAG/{sample}.tsv"
    threads: 8
    envmodules:
        "coverm/0.6.1",
        "samtools/1.9-goolf-1.7.20"
    shell:
        """
        coverm contig \
            --methods count \
            --coupled {input.fq1} {input.fq2} \
            --reference {input.ref} \
            -t {threads} \
            -o {output.tsv} \
            --exclude-supplementary  
        """

rule MAG_create_counts_matrix:
    input:
        fna = expand(config["params"]["output"] + "coverM/MAG/{sample}.tsv",sample=SAMPLE)
    params:
        bin = config["bin"]["scripts"],
        dir = config["params"]["output"]
    output:
        tsv = config["params"]["output"] + "tables/counts.coverM.MAG.tsv"
    envmodules:
        "R/4.2.1"
    shell:
        """
        Rscript {params.bin}/coverM_MAG_create_counts_matrix.R -d {params.dir}
        """

rule coverM_viralMAG:
    input:
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq.gz",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq.gz",
        ref = config["params"]["output"] + "VirSorter2/viral_MAGs.fna"
    params:
        dir = config["params"]["output"] + "VirSorter2/pass2/"
    output:
        tsv = config["params"]["output"] + "coverM/viralMAG/{sample}.tsv"
    threads: 1
    envmodules:
        "coverm/0.6.1",
        "samtools/1.9-goolf-1.7.20"
    shell:
        """
        coverm contig \
            --methods count \
            --coupled {input.fq1} {input.fq2} \
            --reference {input.ref} \
            -t {threads} \
            -o {output.tsv} \
            --exclude-supplementary  
        """

rule viralMAG_create_counts_matrix:
    input:
        fna = expand(config["params"]["output"] + "coverM/viralMAG/{sample}.tsv",sample=SAMPLE)
    params:
        bin = config["bin"]["scripts"],
        dir = config["params"]["output"]
    output:
        tsv = config["params"]["output"] + "tables/counts.coverM.viralMAG.tsv"
    envmodules:
        "R/4.2.1"
    shell:
        """
        Rscript {params.bin}/coverM_viralMAG_create_counts_matrix.R -d {params.dir}
        """

rule BLASTn:
    input:
        fna = config["params"]["output"] + "VirSorter2/pass2/{sample}/final-viral-combined.fa"
    params:
        db = config["ref"]["viral_nt_db"]
    output:
        tsv = config["params"]["output"] + "BLASTn/{sample}.blastn.tsv"
    threads: 8
    envmodules:
        "blast+/2.13.0"
    shell:
        """
        blastn -db {params.db} -query {input.fna} -out {output.tsv} -num_threads {threads} -outfmt 6
        """

rule DRAM_v:
    input:
        fna = config["params"]["output"] + "VirSorter2/pass2/{sample}/for-dramv/final-viral-combined-for-dramv.fa",
        tab = config["params"]["output"] + "VirSorter2/pass2/{sample}/for-dramv/viral-affi-contigs-for-dramv.tab"
    params:
        dir1 = config["params"]["output"] + "dramv/",
        dir2 = config["params"]["output"] + "dramv/annotate/{sample}/",
        dir3 = config["params"]["output"] + "dramv/distill/{sample}/",
        bp = config["params"]["min_MAG_bp_size"]
    output:
        tsv = config["params"]["output"] + "dramv/distill/{sample}/vMAG_stats.tsv"
    threads: 4
    envmodules:
        "dram/1.3.4-Python-3.10.4"
    shell:
        """
        if [ -s {input.fna} ]
        then
            rm -rf {params.dir2}
            DRAM-v.py annotate -i {input.fna} -v {input.tab} -o {params.dir2} --skip_trnascan --threads {threads} --min_contig_size {params.bp}
            if [ $(grep "vogdb_categories" {params.dir2}/annotations.tsv | wc -l) -lt 1 ]; then sed -i "s/amg_flags$/amg_flags\tvogdb_categories/g" {params.dir2}/annotations.tsv; fi

            rm -rf {params.dir3}

            DRAM-v.py distill -i {params.dir2}/annotations.tsv -o {params.dir3}
        else
            mkdir -p {params.dir2}
            mkdir -p {params.dir3}
            touch {output.tsv}
        fi
        """

rule Prodigal:
    input:
        fna = config["params"]["output"] + "VirSorter2/pass2/{sample}/final-viral-combined.fa"
    output:
        genes=config["params"]["output"] + "Prodigal/{sample}.viral_genomes.genes",
        faa=config["params"]["output"] + "Prodigal/{sample}.viral_genomes.faa"
    threads: 8
    envmodules:
        "prodigal/2.6.3"
    shell:
        """
        if [ -s {input.fna} ]
        then
            prodigal -i {input.fna} -o {output.genes} -a {output.faa} -p meta
        else
            touch {output.genes}
            touch {output.faa}
        fi
        """

rule vConTACT2:
    input:
        faa = config["params"]["output"] + "Prodigal/{sample}.viral_genomes.faa"
    params:
        clustalone = config["bin"]["clustalone"],
        dir = config["params"]["output"] + "vConTACT2/{sample}/",
        gene2genome = config["params"]["output"] + "vConTACT2/{sample}/gene2genome.csv"
    output:
        csv = config["params"]["output"] + "vConTACT2/{sample}/viral_cluster_overview.csv"
    threads: 8
    envmodules:
        "vcontact2/0.11.3-Python-3.10.8-Snakemake-7.25.0"
    shell:
        """
        rm -rf {params.dir}
        mkdir {params.dir}
        vcontact2_gene2genome -p {input.faa} -o {params.gene2genome} -s Prodigal-FAA
        vcontact2 --raw-proteins {input.faa} --rel-mode BLASTP --proteins-fp {params.gene2genome} --db ProkaryoticViralRefSeq94-Merged --pcs-mode MCL --vcs-mode ClusterONE --output-dir {params.dir} -t {threads} --c1-bin {params.clustalone}/cluster_one-1.0.jar
        """

rule MetaPhlAn:
    input:
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq.gz",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq.gz"
    params:
        db = config["ref"]["metaphlan_db"]
    output:
        bowtie2out = config["params"]["output"] + "MetaPhlAn/{sample}.bowtie2out.txt",
        txt = config["params"]["output"] + "MetaPhlAn/{sample}.MetaPhlAn.txt"
    threads: 8
    envmodules:
        "metaphlan/4.0.2"
    shell:
        """
        metaphlan \
            {input.fq1},{input.fq2} \
            --bowtie2db {params.db} \
            --bowtie2out {output.bowtie2out} \
            -o {output.txt} \
            --nproc {threads} \
            -t rel_ab_w_read_stats \
            --add_viruses --input_type fastq
        """

rule HUMAnN:
    input:
        fq1 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.1.fq.gz",
        fq2 = config["params"]["output"] + "KrakenTools/{sample}.exclude.s_Homo_sapiens.2.fq.gz"
    params:
        dir = config["params"]["output"] + "HUMAnN/{sample}",
        fq = config["params"]["output"] + "HUMAnN/{sample}/{sample}.fq.gz",
        humann_db = config["ref"]["humann_db"],
        metaphlan_db = config["ref"]["metaphlan_db"]
    output:
        txt = config["params"]["output"] + "HUMAnN/{sample}/test.txt"
    threads: 8
    envmodules:
        "metaphlan/4.0.2",
        "humann/3.6"
    shell:
        """
        cat {input.fq1} {input.fq2} > {params.fq}
                
        humann \
            -i {params.fq} \
            -o {params.dir} \
            --input-format fastq.gz \
            --nucleotide-database {params.humann_db}/chocophlan \
            --protein-database {params.humann_db}/uniref \
            --metaphlan-options "-t rel_ab_w_read_stats --bowtie2db {params.metaphlan_db}" \
            --threads {threads}

        touch {output.txt}
        rm {params.fq}
        """

rule PILER_CR:
    input:
        fa = config["params"]["output"] + "SPAdes/{sample}/contigs.final.fasta"
    output:
        fa = config["params"]["output"] + "PILER-CR/{sample}.direct_repeats.fa",
        txt = config["params"]["output"] + "PILER-CR/{sample}.report.txt"
    envmodules:
        "pilercr/1.0.6"
    shell:
        """
        pilercr -in {input.fa} -out {output.txt} -seq {output.fa} -noinfo
        """

rule PILER_CR_concatenate:
    input:
        fa = expand(config["params"]["output"] + "PILER-CR/{sample}.direct_repeats.fa",sample=SAMPLE)
    output:
        fa = config["params"]["output"] + "PILER-CR/direct_repeats.fa"
    envmodules:
        "pilercr/1.0.6"
    shell:
        """
        cat {input.fa} > {output.fa}
        """

rule MetaCRAST:
    input:
        dr = config["params"]["output"] + "PILER-CR/direct_repeats.fa",
        fa = config["params"]["output"] + "SPAdes/{sample}/contigs.final.fasta"
    params:
        dir = config["params"]["output"] + "MetaCRAST/{sample}"
    output:
        fa = config["params"]["output"] + "MetaCRAST/{sample}/totalSpacers.fa"
    envmodules:
        "MetaCRAST/20200309"
    shell:
        """
        MetaCRAST -p {input.dr} -i {input.fa} -o {params.dir} -d 3
        """

rule MAG_summary:
    input:
        amrfinderplus_output = expand(config["params"]["output"] + "AMRFinderPlus/{sample}.tsv",sample=SAMPLE),
        dramv_output = expand(config["params"]["output"] + "dramv/distill/{sample}/vMAG_stats.tsv",sample=SAMPLE),
        gtdbtk_output = config["params"]["output"] + "GTDB-Tk/classify/gtdbtk.bac120.summary.tsv",
        iphop_output = expand(config["params"]["output"] + "iPHoP/{sample}/Host_prediction_to_genus_m90.csv",sample=SAMPLE),
        kraken2_output = expand(config["params"]["output"] + "Kraken2_MAG/{sample}.output.tsv",sample=SAMPLE),
        vamb_output = config["params"]["output"] + "VAMB/clusters/clusters.tsv",
        vcontact2_output = expand(config["params"]["output"] + "vConTACT2/{sample}/viral_cluster_overview.csv",sample=SAMPLE),
        viral_blastn_output = expand(config["params"]["output"] + "BLASTn/{sample}.blastn.tsv",sample=SAMPLE)
    params:
        bin = config["bin"]["scripts"],
        dir = config["params"]["output"]
    output:
        tsv=config["params"]["output"] + "tables/MAG_summary.tsv"
    envmodules:
        "R/4.2.1"
    shell:
        """
        Rscript {params.bin}/summarize_MAGs.R -d {params.dir}
        """