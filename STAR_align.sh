#!/bin/bash
#SBATCH --ntasks-per-node=12
#SBATCH --time=10:00:00
#SBATCH --mem=40GB
#SBATCH --nodes=1

#Load variable folder
#I use a folder with variables for each sample so it's easier to run
#This can be done in many different ways
. variables

## Create output directory structure
## CHECK PATHS ARE CORRECT !!
fastq_dir="/mainfs/BaralleLab/for_luke/SOT111"
sample_dir="/mainfs/BaralleLab/for_luke/$sampleID"

#Check fastq quality
module load biobuilds/2017.05
fastqc -o /mainfs/BaralleLab/for_Luke/SOT111/ /mainfs/BaralleLab/for_Luke/SOT111/SOT111_R1.merged.fastq.gz /mainfs/BaralleLab/for_Luke/SOT111/SOT111_R2.merged.fastq.gz --threads 12

################################################################################
#Align with STAR
#quantMode GeneCounts > counts number of reads while mapping
#--outReadsUnmapped Fastx > outputs in separate file (fasta) unmapped reads 
#--outFilterType BySJout > reuces the number of "spurious" junctions
#################################################################################

#STAR is available in IRIDS you need to module load and remove the path below
#I have it installed locally
#You will need to create your own STAR index (see documentation)

module load STAR/2.7.10a
STAR --genomeDir /mainfs/BaralleLab/for_Luke/genome_data/star_index/ \
	--readFilesCommand zcat \
	--readFilesIn /mainfs/BaralleLab/for_Luke/SOT111/SOT111_R1.merged.fastq.gz /mainfs/BaralleLab/for_Luke/SOT111/SOT111_R2.merged.fastq.gz \
	--runThreadN 12 \
	--twopassMode Basic \
	--twopass1readsN -1 \
	--outSAMmapqUnique 60 \
	--outFilterType BySJout \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--quantMode GeneCounts \
	--outSAMtype BAM Unsorted

## Clean up STAR intermediates
rm -r _STARgenome
rm -r _STARpass1
	
####################################################
## Samtools (v1.9) sort and index the alignments
####################################################

module load samtools/1.9

## Sort main alignment
samtools sort -@ 12 Aligned.out.bam > "$sampleID".main.sorted.bam
samtools index "$sampleID".main.sorted.bam

## Clean up unsorted bams as they take up a lot of space
rm Aligned.out.bam

