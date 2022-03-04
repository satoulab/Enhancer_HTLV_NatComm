#!/bin/bash

#3_IS_detection.sh

# This script is used to identify virus integration sites

# Software required:
#               - bwa
#               - samtools
#               - picard

# In-house Python scripts:
#               - detect_IS_v2.py
#               - grep_reads.py

# ************************************* SCRIPT BEGINS HERE ************************************* #

# Make analysis directory

mkdir 3_virus_host_junction
cd ./3_virus_host_junction

# Mapping performed using bwa against reference genome hg19_HTLV_LTR
# The reference genome contains all human chromosomes (chr1-22, X, Y, M)
# The virus is separated into 2 different chromosomes - the LTR only and the whole virus excluding LTRs echo "Task started at"

# Mapping to reference genome (hg19+AB513134+LTR)
	echo
	echo "Mapping using bwa"
	echo "Mapping against hg19+HTLV-1+LTR"
		bwa mem -t 4 -M -R "@RG\tID:sample\tSM:sample\tPL:Illumina" /reference_genome/hg19+HTLV-1+LTR/genome.fa ../1_qc/R1_clean.fastq ../1_qc/R2_clean.fastq > hg19+virus+LTR.sam
	echo
	echo "Mapping complete - proceed to extract reads"
	echo

# Reads extraction and selection
	echo "Begin reads extraction and selection"
	echo
		samtools view -Sb hg19+virus+LTR.sam > hg19+virus+LTR.bam
	echo
	echo "Filter and accept only first mapped reads"
		samtools view -bh -F 256 -o hg19+virus+LTR_uniq.bam hg19+virus+LTR.bam
	echo
		samtools sort hg19+virus+LTR_uniq.bam > hg19+virus+LTR_uniq_sort.bam
	echo
     	echo "Duplicates removal"
                picard MarkDuplicates INPUT=./hg19+virus+LTR_uniq_sort.bam OUTPUT=./hg19+virus+LTR_uniq_sort_duprmv.bam METRICS_FILE=marked_dup_metrics.txt TMP_DIR=./tmp REMOVE_DUPLICATES=TRUE MAX_RECORDS_IN_RAM=null ASSUME_SORTED=true

	echo "Convert bam to sam file"
		samtools view -h hg19+virus+LTR_uniq_sort_duprmv.bam > hg19+virus+LTR_uniq_sort_duprmv.sam
	echo
	echo "Filter reads containing HTLV-1"

awk -F "\t" '($3 ~ /AB513134_LTRs/ || $7 ~/AB513134_LTRs/ || $3 ~/AB513134_noLTR/ || $7 ~/AB513134_noLTR/ || $1~/^@/) {print}' hg19+virus+LTR_uniq_sort_duprmv.sam > virus_LTR_duprmv_uniq_sort.sam

# making total_chimera_paired_sort.bam file

awk -F "\t" '($3 ~ /AB513134_LTRs/ && $7 ~/chr/ || $3 ~ /chr/ && $7 ~/AB513134_LTRs/ || $3 ~ /AB513134_noLTR/ && $7 ~/chr/ || $3 ~ /chr/ && $7 ~/AB513134_noLTR/ || $1~/^@/) {print}' virus_LTR_duprmv_uniq_sort.sam > total_chimera_paired.sam
		samtools view -Sb total_chimera_paired.sam > total_chimera_paired.bam

# IS_analysis by using python
        echo "activate python2.7"
               source activate py27_bioinfo
        echo
	echo "perform IS analysis"
        

          python2.7 ../In-house Python scripts/detect_Is.py total_chimera_paired.bam ../txt/hg19+AB513134_LTR+noLTR_chrome_size.txt 1000 20


# extract Total_IS.sam
	  python2.7 ../In-house Python scripts/grep_reads.py total_chimera_paired.bam ./total_chimera_paired.bam.out.txt ./
	echo
	awk -F"\t" ' $1~/^@/ {print}' total_chimera_paired.sam > header
	echo
		cat header IS* > total_IS.sam 
	echo	
		samtools view -Sb total_IS.sam > total_IS.bam
		samtools sort total_IS.bam total_IS_sort
		samtools view -h total_IS_sort.bam > total_IS_sort.sam

# remove intermediate files

	echo "Congratulations! Data processing complete!"
	echo "Task completed on"
		date
	end=$(date +%s)
	duration=$(($end-$begin))
 	echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
	echo
#
# ************ Main Script Ends Here ************ #

