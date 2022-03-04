#!/bin/bash


## This script is used for alignment of reads against reference genome

# Software required:
#               - bwa
#               - samtools
		- picard

# ************************************* SCRIPT BEGINS HERE ************************************* #
#
	echo
	echo "Begin analysis"
	echo "Task started at"
		date
	begin=$(date +%s)

	echo "Proceed to step 2 in 2_map directory"
		cd ./2_map		

# Mapping to reference genome
	echo 
	echo "Mapping using bwa"
	echo "Mapping against hg19+htlv1_AB513134"
		bwa mem -t 4 -Y -L 0 -M -R "@RG\tID:sample\tSM:sample\tPL:Illumina" /reference_genome/hg19+HTLV_AB513134/genome.fa ../1_qc/R1_clean.fastq ../1_qc/R2_clean.fastq > hg19+virus.sam
	echo 
	echo "Mapping complete - proceed to extract reads"
	echo 
# Reads extraction and selection
	echo "Begin reads extraction and selection"
	echo 
		samtools view -Sb hg19+virus.sam > hg19+virus.bam
	echo 
	echo "Filter and accept only first mapped reads"
		samtools view -bh -F 256 -o hg19+virus_uniq.bam hg19+virus.bam
		samtools sort hg19+virus_uniq.bam > hg19+virus_uniq_sort.bam

	echo "Duplicates removal"
                picard MarkDuplicates INPUT=./hg19+virus_uniq_sort.bam OUTPUT=./hg19+virus_uniq_sort_duprmv.bam METRICS_FILE=marked_dup_metrics.txt REMOVE_DUPLICATES=TRUE MAX_RECORDS_IN_RAM=null ASSUME_SORTED=true TMP_DIR=./
	echo "Convert bam file to sam file"
		samtools view -h hg19+virus_uniq_sort_duprmv.bam > hg19+virus_uniq_sort_duprmv.sam
	echo 
	echo "Filter reads which map to HTLV-1"
		awk -F"\t" '($3 ~ /AB513134/ || $7 ~/AB513134/ || $1~/^@/) {print}' hg19+virus_uniq_sort_duprmv.sam > virus_uniq_map_duprmv.sam
	echo
		awk -F"\t" '($5>20 || $1~/^@/) {print}' virus_uniq_map_duprmv.sam > virus_uniq_map_duprmv_MAPQ20.sam

	echo 
	echo "Extract reads with soft-clipping"
        	awk -F"\t" '($6 ~/S/ || $1~/^@/) {print}' virus_uniq_map_duprmv.sam > total_softclipping.sam
	echo "sam-to-bam"
		samtools view -Sb virus_uniq_map_duprmv.sam > virus_uniq_map_duprmv.bam
        echo
# gzip fastq file 
		gzip ../R1*
		gzip ../R2*
		gzip ../I*	

# Count number of reads 
	echo "Analysis complete - displaying read counts"

	echo "Mapped reads"
		samtools view -F 0x4 ../2_map/hg19+virus_uniq_sort_duprmv.bam | cut -f 1 | sort | uniq | wc -l

	echo "Reads mapped to HTLV-1"
		samtools view -F 0x4 ../2_map/virus_uniq_map_duprmv.bam | cut -f 1 | sort | uniq | wc -l

	echo "Soft-clipping Reads mapped to total viral reads before cleaning"
                samtools view -S -F 0x4 ../2_map/total_softclipping.sam | cut -f 1 | sort | uniq | wc -l

	 "Congratulations! Data processing complete!"
	echo "Task completed on"
		date
	end=$(date +%s)
	duration=$(($end-$begin))
	echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
	echo
#
# ************ Main Script Ends Here ************ #
