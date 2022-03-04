#!/bin/bash

# 1_Adaptor_Remove_Index_QC.sh

# This script first removes reads with low quality index reads followed by removing the sequencing
# adaptor sequence from the reads
# Finally, it performs another quality check on the reads using FASTQC

# Software required:
#               - cutadapt

# In-house Perl scripts:
#               - IndexQuality_I8_20.pl
#               - qcleaner_renew_v3.1.pl

# ************************************* SCRIPT BEGINS HERE ************************************* #
# ************ Main Script Begins Here ************ #
        echo
	echo
	echo "Begin analysis"
	echo "Task started at"
		date
	begin=$(date +%s)

# Prepare fastq files
	echo "Unzip"
		gunzip *gz
	
	echo "convert fastq file name"
		mv  *R1* R1.fastq
		mv  *R2* R2.fastq
		mv  *I1* I1.fastq

# Make analysis directories
	echo
	echo "Create analysis directories"
		mkdir 1_qc
		mkdir 2_map

	echo "Directories 1_qc & 2_map created"
	echo "Begin step 1 in 1_qc directory"
		cd 1_qc

# Remove reads with low quality of Index (iQC20)
	echo
	echo " Remove reads with low Index Quality (iQC20)"
	echo
		perl ../In-house Perl scripts/IndexQuality_I8_20.pl ../R1.fastq ../R2.fastq ../I1.fastq ./R1_20.fastq ./R2_20.fastq ./I1_20.fastq
	
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
		../In-house Perl scripts/qcleaner_renew_v3.1.pl --i1 R1_step1.fastq --i2 R2_step1.fastq --o1 R1_clean.fastq --o2 R2_clean.fastq --log qclog.txt --trim skip
	
# Cleaning intermediated files
	rm R1_step1.fastq
	rm R2_step1.fastq
	rm R1_20.fastq
	rm R2_20.fastq
	rm I1_20.fastq 
        echo
        echo
	 "Congratulations! Data processing complete!"
	echo "Task completed on"
		date
	end=$(date +%s)
	duration=$(($end-$begin))
	echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
	echo
#
# ************ Main Script Ends Here *******#
