#!/bin/bash

####################################################
#                                                  # 
#  	       QC_MAP_hg19+HTLV_SE.sh 		   #
#                       		           #  
#						   #
#						   #
####################################################

# ************ Main Script Begins Here ************ #

	echo
	echo "Begin analysis"
	echo "Task started at"
                date
        begin=$(date +%s)

# Make analysis directories
	echo
	echo "Create analysis directories"
		mkdir 1_qc
		mkdir 2_map
	echo "Combine fastq file"
		cat *.fastq.gz > R1.fastq.gz
		gunzip R1.fastq.gz
		rm *R1_001.fastq.gz
	echo
        echo
                cd 1_qc
        echo
# Remove adaptor sequence
	echo
	echo "Remove adaptor from read 1"
	echo
		cutadapt -m 1 -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 5 ../R1.fastq -o R1_combined_step1.fastq
	echo

# Read QC
	echo "Reads quality check"
	echo
		./prinseq-lite.pl -min_len 10 -trim_qual_right 20 -min_qual_mean 20 -fastq ./R1_combined_step1.fastq -out_good clean_R
	
# Cleaning intermediated files
	rm R1_combined_step1.fastq

	echo "Proceed to step 2 in 2_map directory"
		cd ../2_map		

# Mapping to reference genome
	echo 
	echo "Mapping using bwa"
	echo "Mapping against hg19+htlv1_AB513134"
		bwa mem -t 4 -Y -L 0 -M -R "@RG\tID:sample\tSM:sample\tPL:Illumina" /reference_genome/BWA/HTLV/hg19+HTLV_AB513134/genome.fa ../1_qc/clean_R.fastq > hg19+virus.sam
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
		samtools sort hg19+virus_uniq.bam -o hg19+virus_uniq_sort.bam

	echo "Duplicates removal"
                java -jar /programs/picard.jar MarkDuplicates INPUT=./hg19+virus_uniq_sort.bam OUTPUT=./hg19+virus_uniq_sort_duprmv.bam METRICS_FILE=marked_dup_metrics.txt REMOVE_DUPLICATES=TRUE ASSUME_SORTED=true
	echo "Convert bam file to sam file"
		samtools view -h hg19+virus_uniq_sort_duprmv.bam > hg19+virus_uniq_sort_duprmv.sam
		samtools index hg19+virus_uniq_sort_duprmv.bam 
	echo 
	echo "Filter reads which map to HTLV-1"
		awk -F"\t" '($3 ~ /AB513134/ || $7 ~/AB513134/ || $1~/^@/) {print}' hg19+virus_uniq_sort_duprmv.sam > virus_uniq_map_duprmv.sam
	echo
		awk -F"\t" '($5>20 || $1~/^@/) {print}' virus_uniq_map_duprmv.sam > virus_uniq_map_duprmv_MAPQ20.sam

	echo 
	echo "sam-to-bam"
		samtools view -Sb virus_uniq_map_duprmv.sam > virus_uniq_map_duprmv.bam
		samtools view -Sb virus_uniq_map_duprmv_MAPQ20.sam > virus_uniq_map_duprmv_MAPQ20.bam
		samtools index virus_uniq_map_duprmv.bam
		samtools index virus_uniq_map_duprmv_MAPQ20.bam

# Cleaning intermediated files
		rm hg19+virus.sam
		rm hg19+virus.bam
		rm hg19+virus_uniq_sort.bam
		rm hg19+virus_uniq.bam
		rm hg19+virus_uniq_sort_duprmv.sam
		rm virus_uniq_map_duprmv_MAPQ20.sam	

# Count number of reads 
	echo "Analysis complete - displaying read counts"

	echo "Mapped reads"
		samtools view -F 0x4 ../2_map/hg19+virus_uniq_sort_duprmv.bam | cut -f 1 | sort | uniq | wc -l

	echo "Reads mapped to HTLV-1"
		samtools view -F 0x4 ../2_map/virus_uniq_map_duprmv.bam | cut -f 1 | sort | uniq | wc -l

# Compress fastq file
	echo "compression"
		gzip ../R1.fastq
		gzip ../1_qc/clean_R.fastq

	 "Congratulations! Data processing complete!"
	echo "Task completed on"
		date
	end=$(date +%s)
	duration=$(($end-$begin))
	echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
	echo
#
# ************ Main Script Ends Here ************ #
