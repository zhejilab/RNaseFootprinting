#!/bin/bash


## Sample description: SRX8494468, GSM4594590, HeLa_RFP_gCtrl1; Homo sapiens; RNA-Seq
SAMPLE="SRR11950060"


## Step 1: download reference files

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


## Step 2: make indices for Bowtie2 and TopHat

module purge all
module load bowtie2/2.2.6
module load tophat/2.1.0

mkdir ./annotation/rRNA
bowtie2-build -f ./annotation/rRNA.fa ./annotation/rRNA/rRNA
mv ./annotation/rRNA.fa ./annotation/rRNA/rRNA.fa

mkdir ./annotation/genome
bowtie2-build -f ./annotation/genome.GRCh38.fa ./annotation/genome/genome
mv ./annotation/genome.GRCh38.fa ./annotation/genome/genome.fa

mkdir ./annotation/transcriptome
tophat -G ./annotation/gencode.v42.gtf --transcriptome-index=./annotation/transcriptome/known ./annotation/genome/genome
mv ./annotation/gencode.v42.gtf ./annotation/transcriptome/known.gtf


## Step 3: download FASTQ file from Gene Expression Omnibus repository

mkdir ./data
module purge all
module load sratoolkit/3.0.0
fasterq-dump --outdir ./data/ ${SAMPLE}

gzip ./data/${SAMPLE}_1.fastq
gzip ./data/${SAMPLE}_2.fastq


## Step 4: trim adapters from input reads

cutadapt -a AAAAAAAAAA -e 0.1 -u 7 -m 18 -M 35 \
 -o ./data/${SAMPLE}_1.trimmed.fastq.gz \
 ./data/${SAMPLE}_1.fastq.gz > ./logs/${SAMPLE}_cutadapt.log


## Step 5: align trimmed reads to rRNA reference sequences using Bowtie2

mkdir ./align_rRNA
module purge all
module load bowtie2/2.2.6

bowtie2 -p 8 \
 -x ./annotation/rRNA/rRNA \
 -U ./data/${SAMPLE}_1.trimmed.fastq.gz \
 --un ./align_rRNA/unmapped.fastq \
 -S ./align_rRNA/accepted_hits.sam 2> ./logs/${SAMPLE}_align_rRNA.log


## Step 6: align rRNA-unmapped reads to the transcriptome using TopHat

mkdir ./align_transcriptome
module purge all
module load python/anaconda
module load bowtie2/2.2.6
module load tophat/2.1.0

tophat -p 8 --keep-tmp \
 --transcriptome-index=./annotation/transcriptome/known \
 -o ./align_transcriptome/ \
 ./annotation/genome/genome \
 ./align_rRNA/unmapped.fastq &> ./logs/${SAMPLE}_align_transcriptome.log


## Step 7: create BigWig tracks showing the genome coverage of uniquely-mapped reads 

module purge all
module load deeptools/3.1.1
module load samtools/1.6

samtools view -H ./align_transcriptome/accepted_hits.bam > ./align_transcriptome/unique_hits.sam
samtools view ./align_transcriptome/accepted_hits.bam | grep -w "NH:i:1" >> ./align_transcriptome/unique_hits.sam
samtools view -S -b ./align_transcriptome/unique_hits.sam | samtools sort -o ./align_transcriptome/unique_hits.bam
samtools index ./align_transcriptome/unique_hits.bam
bamCoverage -p 8 -b ./align_transcriptome/unique_hits.bam -o ./align_transcriptome/unique_hits.str1.bw --filterRNAstrand reverse --normalizeUsing CPM --binSize 1
bamCoverage -p 8 -b ./align_transcriptome/unique_hits.bam -o ./align_transcriptome/unique_hits.str2.bw --filterRNAstrand forward --normalizeUsing CPM --binSize 1


## Step 8: quantify uniquely-mapped reads in known transcripts using HT-Seq

htseq-count \
 --format="bam" \
 --stranded="yes" \
 --type="CDS" \
 --mode="union" \
 --nonunique="all" \
 ./align_transcriptome/unique_hits.bam \
 ./annotation/transcriptome/known.gff > ./align_transcriptome/unique_hits.counts.txt


## Step 9: annotate the reference genome annotation and primary assembly to identify candidate ORFs

mkdir -p ./riborf/annotate
perl ./software/RibORF.2.0/ORFannotate.pl -g ./annotation/genome/genome.fa -t ./annotation/refseq.genePred -o ./riborf/annotate


## Step 10: quality control analyses to plot the distribution of reads surrounding the start and stop codons by read length

mkdir -p ./riborf/readdist
module purge all
module load R/3.6.3

perl ./software/RibORF.2.0/readDist.pl \
 -f ./align_transcriptome/unique_hits.sam \
 -g ./annotation/refseq.genePred \
 -o ./riborf/readdist \
 -d 18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,

Rscript --vanilla software/RibORF.2.0/combine_plots.R ./riborf/readdist 18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,


## Step 11: correct the location of reads to ribosomal A-sites based on their read length and the provided offset parameter file that needs to be manually curated

mkdir -p ./riborf/corrected
module purge all
module load samtools/1.6
module load R/3.6.3

perl ./software/RibORF.2.0/offsetCorrect.pl \
 -r ./align_transcriptome/unique_hits.sam \
 -p ./riborf/offset.correction.parameters.txt \
 -o ./riborf/corrected/corrected_hits.sam

perl ./software/RibORF.2.0/readDist.pl \
 -f ./riborf/corrected/corrected_hits.sam \
 -g ./annotation/refseq.genePred \
 -o ./riborf/corrected \
 -d 1

Rscript --vanilla ./software/RibORF.2.0/combine_plots.R ./riborf/corrected 1

samtools view -S -b ./riborf/corrected/corrected_hits.sam | samtools sort -o ./riborf/corrected/corrected_hits.bam
samtools index ./riborf/corrected/corrected_hits.bam
bamCoverage -p 8 -b ./riborf/corrected/corrected_hits.bam -o ./riborf/corrected/corrected_hits.str1.bw --filterRNAstrand reverse --normalizeUsing CPM --binSize 1
bamCoverage -p 8 -b ./riborf/corrected/corrected_hits.bam -o ./riborf/corrected/corrected_hits.str2.bw --filterRNAstrand forward --normalizeUsing CPM --binSize 1


