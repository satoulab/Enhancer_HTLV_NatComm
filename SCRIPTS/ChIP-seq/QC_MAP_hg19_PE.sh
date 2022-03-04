#!/bin/bash

####################################################
#                                                  # 
#    	         QC_MAP_hg19_PE.sh      	   #
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
	echo
# Decompress and combind
        echo
               gunzip *.fastq.gz
	echo "convert fastq file name"
		mv  *R1* R1.fastq
		mv  *R2* R2.fastq
        echo
                cd 1_qc
        echo

# Remove adaptor sequence
	echo  
	echo "Remove adaptor from read 1"
	echo 
		cutadapt -m 1 -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 5 ./R1_20.fastq -o R1_step1.fastq
	echo 
	echo "Remove adaptor from read 2"
	echo 
		cutadapt -m 1 -b GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 ./R2_20.fastq -o R2_step1.fastq


# Read QC
	echo "Reads quality check"
	echo
		./prinseq-lite.pl -min_len 10 -trim_qual_right 20 -min_qual_mean 20 -fastq ./R1_step1.fastq -fastq2 ./R2_step1.fastq -out_good clean_R
	
# Cleaning intermediated files
	rm R1_step1.fastq
	rm R2_step1.fastq
	rm R1_20.fastq
	rm R2_20.fastq

	echo "Proceed to step 2 in 2_map directory"
		cd ../2_map		

# Mapping to reference genome
	echo 
	echo "Mapping using bwa"
	echo "Mapping against hg19"
		bwa mem -t 4 -Y -L 0 -M -R "@RG\tID:sample\tSM:sample\tPL:Illumina" /reference_genome/hg19/genome.fa ../1_qc/clean_R_1.fastq ../1_qc/clean_R_2.fastq > hg19.sam
	echo 
	echo "Mapping complete - proceed to extract reads"
	echo 
# Reads extraction and selection
	echo "Begin reads extraction and selection"
	echo 
		samtools view -Sb hg19.sam > hg19.bam
	echo 
	echo "Filter and accept only first mapped reads"
		samtools view -bh -F 256 -o hg19_uniq.bam hg19.bam
		samtools sort hg19_uniq.bam -o hg19_uniq_sort.bam

	echo "Duplicates removal"
                java -jar /programs/picard.jar MarkDuplicates INPUT=./hg19_uniq_sort.bam OUTPUT=./hg19_uniq_sort_duprmv.bam METRICS_FILE=marked_dup_metrics.txt REMOVE_DUPLICATES=TRUE ASSUME_SORTED=true
	echo "Generate index"
		samtools index hg19_uniq_sort_duprmv.bam 
	echo 
# Cleaning intermediated files
		rm hg19.sam
		rm hg19.bam
		rm hg19_uniq_sort.bam
		rm hg19_uniq.bam	

# Count number of reads 
	echo "Analysis complete - displaying read counts"

	echo "Mapped reads"
		samtools view -F 0x4 ../2_map/hg19_uniq_sort_duprmv.bam | cut -f 1 | sort | uniq | wc -l
	echo
	echo	 "Congratulations! Data processing complete!"
	echo "Task completed on"
		date
	end=$(date +%s)
	duration=$(($end-$begin))
	echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
	echo
#
# ************ Main Script Ends Here ************ #
