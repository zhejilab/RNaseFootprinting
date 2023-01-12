# Rfoot-seq Protocol and Analysis

This repository contains the code for analyzing Rfoot-seq data, a platform for rapid ribosome profiling in low-input samples. Rfoot-seq can be used to map stable RNA-protein complexes using RNase digestion without complex crosslinking. This allows for the study of simultaneous cytosolic and mitochondrial translation via the quantification of ribosome-protected RNA fragments in coding regions. Non-ribosomal RNA-protein complexes associated with non-coding RNAs are also detected and can be used to characterize functional non-coding RNA domains.

Contact *zhe.ji (at) northwestern.edu* for any questions about the protocol or analysis.


## Analysis

### Rfoot-seq analysis requires the following packages to be installed
- sratoolkit/3.0.0
- python/anaconda
- cutadapt
- bowtie2/2.2.6
- tophat/2.1.0
- deeptools/3.1.1
- samtools/1.6
- htseq
- R/3.6.3
- RibORF.2.0 (installed into the `./software/` directory, available from the [Zhe Ji lab github](https://github.com/zhejilab/RibORF))

### Rfoot-seq analysis requires the following reference files

Prior to running the example script below, these reference files should be placed in the `./annotation/` directory. Note that the file names will need to match those in the list or changed within the example script file.

- `genome.GRCh38.fa` = FASTA file containing genomic sequences split by chromosome from UCSC (downloaded during analysis)
- `gencode.v42.gtf` = GTF file containing genomic annotation from GENCODE (downloaded during analysis)
- `refseq.genePred` = GenePred file containing genomic annotation from RefSeq (downloaded during analysis)
- `rRNA.fa` = FASTA file containing rRNA sequences (available at XX)


### Example analysis for sample SRR11950060

The `example` directory contains a representative end-to-end analysis for sample SRR11950060 with all necessary command line prompts. The sample SRR number can be replaced and run for any sample of interest. Note that the Bowtie2 and TopHat index creation step only needs to be run once.

