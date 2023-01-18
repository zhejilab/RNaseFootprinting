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

Note that these packages will need to be locatable in your PATH variable, or loaded to your workspace if using a computing cluster, before running the commands below.


### Rfoot-seq analysis requires the following reference files

Prior to running the example script below, these reference files should be placed in the `./annotation/` directory. Note that the file names will need to match those in the list or changed within the example script file.

- `genome.GRCh38.fa` = FASTA file containing genomic sequences split by chromosome from UCSC (downloaded during analysis)
- `gencode.v42.gtf` = GTF file containing genomic annotation from GENCODE (downloaded during analysis)
- `refseq.genePred` = GenePred file containing genomic annotation from RefSeq (downloaded during analysis)
- `rRNA.fa` = FASTA file containing rRNA sequences (available in `example/annotation/` directory)


### Example analysis for sample SRR11950060

The `example` directory contains a representative end-to-end analysis for sample SRR11950060 with all necessary command line prompts. The sample SRR number can be replaced and run for any sample of interest. Note that the Bowtie2 and TopHat index creation step only needs to be run once.


#### Read mapping

##### Step 1: Prepare reference files

Download reference files from the UCSC Genome Browser and GENCODE databases. The "refFlat.txt.gz" file contains the gene structure information for each transcript and is processed into genePred format to meet input requirements. All reference files are moved to the annotation directory. These files are needed to create alignment indices for Bowtie2 and TopHat.

```sh
mkdir ./annotation

wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz"
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa ./annotation/genome.GRCh38.fa

wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz"
gunzip gencode.v42.annotation.gtf.gz
mv gencode.v42.annotation.gtf ./annotation/gencode.v42.gtf

wget "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz"
gunzip refFlat.txt.gz
cut refFlat.txt -f2-11 > ./annotation/refseq.genePred
rm refFlat.txt
```

##### Step 2: Make alignment indices for Bowtie2 and TopHat

```sh
mkdir ./annotation/rRNA
bowtie2-build -f ./annotation/rRNA.fa ./annotation/rRNA/rRNA
mv ./annotation/rRNA.fa ./annotation/rRNA/rRNA.fa

mkdir ./annotation/genome
bowtie2-build -f ./annotation/genome.GRCh38.fa ./annotation/genome/genome
mv ./annotation/genome.GRCh38.fa ./annotation/genome/genome.fa

mkdir ./annotation/transcriptome
tophat -G ./annotation/gencode.v42.gtf --transcriptome-index=./annotation/transcriptome/known ./annotation/genome/genome
mv ./annotation/gencode.v42.gtf ./annotation/transcriptome/known.gtf
```

##### Step 3: [Optional] Download FASTQ file of interest from the NIH Gene Expression Omnibus (GEO) repository.

For reproducibility, here we analyze the data for sample "SRR11950060" as an example. This step is optional and the analysis can be running using your own FASTQ file if desired. Just substitute all references to "SRR11950060" with the name of your FASTQ file before running.

```sh
mkdir ./data

fasterq-dump --outdir ./data/ SRR11950060
gzip ./data/SRR11950060_1.fastq
gzip ./data/SRR11950060_2.fastq
```

##### Step 4: Trim read adapters

As the RNase footprints are generally short, the R1 sequencing reads can cover the whole footprint. We used only the R1 reads for the downstream analyses. The first 7 nucleotides are trimmed from the beginning of each read. The polyA adapter is removed from the end of reads, allowing for 1 mismatch within the adapter sequence. Reads with a length between 18-35 nt are retained for downstream analysis to facilitate quick and accurate alignment.

This step takes the downloaded R1 reads as input, performs the read trimming, and outputs another gzipped-FASTQ file containing only the reads that were between 18-35 nt long after removing the 5' and 3' adapters with the specified pattern and error rates.

```sh
cutadapt -a AAAAAAAA -e 0.2 -u 7 -m 18 -M 35 \
 -o ./data/SRR11950060_1.trimmed.fastq.gz \
 ./data/SRR11950060_1.fastq.gz > ./logs/SRR11950060_cutadapt.log
```

##### Step 5: Align trimmed reads to rRNA

Trimmed input reads are aligned to rRNA reference sequences. The reads unmapped to rRNA sequences are used for further analyses. Here we instruct Bowtie2 to output unmapped reads in FASTQ format so that we can use those as the input for the next step, effectively removing any reads that map to rRNA sequences and would reduce our usable read fraction. The other Bowtie2 outputs can be disregarded.

```sh
mkdir ./align_rRNA

bowtie2 -p 8 \
 -x ./annotation/rRNA/rRNA \
 -U ./data/SRR11950060_1.trimmed.fastq.gz \
 --un ./align_rRNA/unmapped.fastq \
 -S ./align_rRNA/accepted_hits.sam 2> ./logs/SRR11950060_align_rRNA.log
```

##### Step 6: Align rRNA-unmapped reads to transcriptome and then genome

We take the rRNA-unmapped reads from the previous step and align them to the reference transcriptome and then genome using TopHat. To prevent an error during alignment, the output directory must be specified with a trailing forward slash.

```sh
mkdir ./align_transcriptome

tophat -p 8 \
 --transcriptome-index=./annotation/transcriptome/known \
 -o ./align_transcriptome/ \
 ./annotation/genome/genome \
 ./align_rRNA/unmapped.fastq &> ./logs/SRR11950060_align_transcriptome.log
```

##### Step 7: Extract uniquely mapped reads and quality check results

We then select only the uniquely mapped reads and generate coverage tracks for quality checking. Coverage tracks in BigWig format for the forward and reverse strands can be viewed in the Integrative Genomics Viewer (IGV).

Extract uniquely mapped reads:

```sh
samtools view -H ./align_transcriptome/accepted_hits.bam > ./align_transcriptome/unique_hits.sam
samtools view ./align_transcriptome/accepted_hits.bam | grep -w "NH:i:1" >> ./align_transcriptome/unique_hits.sam
samtools view -S -b ./align_transcriptome/unique_hits.sam | samtools sort -o ./align_transcriptome/unique_hits.bam
samtools index ./align_transcriptome/unique_hits.bam
```

Generate coverage tracks:

```sh
bamCoverage -p 8 \
 -b ./align_transcriptome/unique_hits.bam \
 -o ./align_transcriptome/unique_hits.str1.bw \
 --filterRNAstrand reverse \
 --normalizeUsing CPM \
 --binSize 1

bamCoverage -p 8 \
 -b ./align_transcriptome/unique_hits.bam \
 -o ./align_transcriptome/unique_hits.str2.bw \
 --filterRNAstrand forward \
 --normalizeUsing CPM \
 --binSize 1
```

#### Data quality control using RibORF

##### Step 8: Examine the read fragment sizes and plot the distribution of 5' read end locations surrounding the start and stop codons grouped by read length.

High quality data should show strong read enrichment in coding regions compared to UTRs. The coding region reads should be enriched fro lengths of ~21 nt and ~29 nt with 3-nt periodicity. Here we use the readDist utility from RibORF to inspect the periodicity and enrichment of different reads by fragment length. The outputs from this step should be used to curate an "offset.correction.parameters.txt" file before running Step 10.

```sh
mkdir -p ./riborf/readdist

perl ./software/RibORF.2.0/readDist.pl \
 -f ./align_transcriptome/unique_hits.sam \
 -g ./annotation/refseq.genePred \
 -o ./riborf/readdist \
 -d 18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,

Rscript --vanilla ./software/RibORF.2.0/combine_plots.R ./riborf/readdist 18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,
```

#### Quantifying ribosome occupancy in protein-coding regions

##### Step 9: Quantify ribosome occupancy in coding regions using HTSeq

We quantify only reads that overlap CDS regions according to the annotation file (--type "CDS"). Reads that overlap a protion of a gene but are not fully contained within that gene are also counted (--mode "union"). Reads that overlap multiple transcript isoforms are counted toward the gene expression (--nonunique "all"). The data is processed in a strand-specific manner (--stranded "yes").

```sh
htseq-count \
 --format="bam" \
 --stranded="yes" \
 --type="CDS" \
 --mode="union" \
 --nonunique="all" \
 ./align_transcriptome/unique_hits.bam \
 ./annotation/transcriptome/known.gff > ./align_transcriptome/unique_hits.counts.txt
```

#### Read distribution features separate ribosomal vs non-ribosomal footprints

##### Step 10: Correct the location of reads based on their read length and the provided offset parameters

The offset parameters will need to be manually curated beefore this step. This will assign the read mapping locations to ribosomal A-sites, correcting for the offset distance between the 5' end of reads and ribosomal A-sites. These correction parameters are inferred from the RibORF readDist plots. Details are described in a [previously published protocol](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpmb.67). Users can load corrected reads (BigWig files) to IGV. The translated ORFs show 3-nt periodicity. And the non-ribosomal footprints show a localized distribution without periodicity.

Correct for read-ribosomal offset distances:

```sh
mkdir -p ./riborf/corrected

perl ./software/RibORF.2.0/offsetCorrect.pl \
 -r ./align_transcriptome/unique_hits.sam \
 -p ./riborf/offset.correction.parameters.txt \
 -o ./riborf/corrected/corrected_hits.sam
```

Confirm correction using readDist visualization:

```sh
perl ./software/RibORF.2.0/readDist.pl \
 -f ./riborf/corrected/corrected_hits.sam \
 -g ./annotation/refseq.genePred \
 -o ./riborf/corrected \
 -d 1

Rscript --vanilla ./software/RibORF.2.0/combine_plots.R ./riborf/corrected 1
```

Visualize corrected reads in IGV:

```sh
samtools view -S -b ./riborf/corrected/corrected_hits.sam | samtools sort -o ./riborf/corrected/corrected_hits.bam
samtools index ./riborf/corrected/corrected_hits.bam

bamCoverage -p 8 \
 -b ./riborf/corrected/corrected_hits.bam \
 -o ./riborf/corrected/corrected_hits.str1.bw \
 --filterRNAstrand reverse \
 --normalizeUsing CPM \
 --binSize 1

bamCoverage -p 8 \
 -b ./riborf/corrected/corrected_hits.bam \
 -o ./riborf/corrected/corrected_hits.str2.bw \
 --filterRNAstrand forward \
 --normalizeUsing CPM \
 --binSize 1
```
