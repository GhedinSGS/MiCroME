# SGS Lab Metagenomics Analysis Pipeline

![image](https://github.com/GhedinSGS/SGSlab_metagenomics/assets/29619358/12df628d-ddd5-47b2-bcef-fab0a922fbd0)

## Overview
The SGS lab metagenomics analysis pipeline is a comprehensive Snakemake pipeline for short paierd-end read metagenomics and metatranscriptomics datasets. The pipeline includes the following steps:

1. **Upstream processing**  
  a. Adaptor trimming, removal of low quality positions, and PCR duplicates using *fastp*  
  b. Mapping-based depletion of host reads using *Bowtie2*  
  c. Additional depletion of human reads using *Kraken2* and *KrakenTools*  
2. **Read-level taxonomic identification and quantification**   
  a. Taxonomic identification and quantification of reads using *Kraken2* and *Bracken* 
3. **Read-level functional annotation and quantification**  
  a. Functional annotation of reads at the KEGG orthology pathway, module, and ortholog level using *FMAP*  
  b. Quantification of antimicrobial resistance genes present in the CARD and MEGARes databases using *CoverM*  
  c. Assembly of metagenome-assembled genomes (MAGs) from processed paired-end reads using *metaSPAdes*  
4. **General in-depth analysis of MAGs**  
  a. Taxonomic identification of MAGs using *Kraken2*  
  b. Annotation of antimicrobial resistance genes in MAGs using *AMRFinderPlus*  
  c. Bacterial-specific taxonomic identification of MAGs using *VAMB* and *GTDB-tk*  
  d. Quantification of MAGs using *coverM*  
  e. Identification of CRISPR repeats and spacers using *PILER-CR* and *MetaCRAST*  
5. **Viral-specific analysis of MAGs**  
  a. Identification of both DNA and RNA viral genomes in MAGs using *VirSorter2*  
  b. Assesss quality and remove host regions for viral MAGs using *CheckV*  
  c. Quantification of viral MAGs using *coverM*  
  d. Annotation of protein-coding genes in viralMAGs using *Prodigal*  
  e. Annotation of viral-specific antimicrobial genes using *DRAM-v*
  f. Identification of bacterial hosts for prophage sequences using *iPHoP*  
  f. Taxonomic identification of viral MAGs using *BLASTn* against viral RefSeq and *vConTACT2*  

## Usage

The pipeline runs using Snakemake, using inputs from the config.yaml. 

```
params:
  reads_dir: "~/reads/"
  output: "~/pipeline_output/"
  min_MAG_bp_size: "500"

bin:
  clustalone: /sysapps/cluster/software/Anaconda3/2022.05/envs/vContact2-0.11.3/bin
  fmap: ~/db/FMAP
  krakentools: ~/packages/KrakenTools_v1.2
  scripts: ~/pipelines/SGSlab_metagenomics/scripts
ref:
  host_ref: ~/reference_files/mouse/GCF_000001635.27_GRCm39_genomic.fna
  card_db: ~/databases/SGSlab_metagenomics/references/card/nucleotide_fasta_protein_variant_model.fasta
  humann_db: ~/databases/SGSlab_metagenomics/references/humann
  kraken2_db: ~/databases/kraken2/20220616_hbvae
  iphop_db: ~/databases/iphop/Sept_2021_pub_rw
  megares_db: ~/databases/megares_full_database_v2.00.fasta
  metaphlan_db: ~/databases/SGSlab_metagenomics/references/metaphlan/mpa_vOct22_CHOCOPhlAnSGB_202212
  viral_nt_db: ~/databases/blast_db/viral.1.1.genomic.fna
```

Under params, the user needs to supply (1) the input directory containing the raw reads for the metagenomics/metatranscriptomics project, (2) an output directory, and (3) the minimum MAG size to be used for downstream MAG-based analyises (recommended 1000 bp for metagenomics studies and 500 bp for metatranscriptomics studies).  

For the initial set up of dependencies, the user must provide directories for the included scripts and the ClustalOne, FMAP, and KrakenTools executables. Additionally, reference genomes and databases for the host organism along with the CARD HUMAnN, Kraken2, iPHoP, MEGARes, MetaPhlAn, and viral RefSeq databases must be provided.

## Output

The pipeline outputs several files including:

1. Counts tables for:  
  a. Read level taxonomic classifications from *Bracken*  
  b. Read level functional annotations from *FMAP*  
  c. Read level antimicrobial gene quantification using CoverM against the *CARD* and *MEGARes* databases  
  d. Assembled MAGs from *metaSPAdes*  
  e. Downstream processed viral MAGs from the viral analysis portion of the pipeline  
2. Annotation tables for:  
  a. CRISPR repeat and spacers identified in MAGs from *PILER-CR* and *MetaCRAST*  
  b. Taxonomic information for assembled MAGs from:  
    -  *VAMB* and *GTDB-Tk*  
    -  *iPHoP*  
    -  *BLASTn* with viral RefSeq  
    -  *vConTACT2* 

## Credits

(add credits including personal citation)  

## Citations

The pipeline uses the following tools:  

(add citations)

