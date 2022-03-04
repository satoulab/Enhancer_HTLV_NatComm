#!/bin/bash

####################################################
#                                                  # 
#            STAR_hg19_HTLV_rmdup_PE.sh            #
#                                                  #  
#                                                  #
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
        echo
               cat *_R1_001.fastq > R1.fastq
               cat *_R2_001.fastq > R2.fastq   
        echo
                cd 1_qc
        echo
# Remove adaptor sequence
	echo
	echo "Remove adaptor from read 1"
	echo
		cutadapt -m 1 -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 5 ../R1.fastq -o R1_combined_step1.fastq
	echo
	echo "Remove adaptor from read 2"
	echo
		cutadapt -m 1 -b GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 ../R2.fastq -o R2_combined_step1.fastq

# Read QC
	echo "Reads quality check"
	echo
		./prinseq-lite.pl -min_len 10 -trim_qual_right 20 -min_qual_mean 20 -fastq ./R1_combined_step1.fastq -fastq2 ./R2_combined_step1.fastq -out_good clean_R
	echo
	echo "Step 1 completed"
	echo
	echo "Proceed to step 2 in 2_map directory"
		cd ../2_map

# Mapping to reference genome
	echo
	echo "Mapping using STAR"
	echo "Mapping against hg19+HTLV-1_AB513134"
		STAR --runThreadN 4 --genomeDir /reference_genome/hg19_AB513134/ --readFilesIn ../1_qc/clean_R_1.fastq ../1_qc/clean_R_2.fastq --chimSegmentMin 20 --outSAMtype BAM SortedByCoordinate
	echo
	echo "Mapping complete - proceed to extract reads"
	echo

# Reads extraction and selection
    echo "Begin reads extraction and selection"
    echo
    echo "Mark duplicate"
        STAR --runMode inputAlignmentsFromBAM --inputBAMfile ./Aligned.sortedByCoord.out.bam --bamRemoveDuplicatesType UniqueIdentical
    echo
    echo "Remove duplicates"
        samtools view -bh -F 0x400 -o Processed_rmdup.bam Processed.out.bam
    echo
    echo "Filter and accept only first mapped reads"
        samtools view -bh -F 256 -o accepted_hits_sort_uniq.bam Processed_rmdup.bam
    echo
    echo "Convert bam file to sam file"
        samtools view -h accepted_hits_sort_uniq.bam > accepted_hits_sort_uniq.sam
    echo
    echo "Count Chimeric.out.junction"
         cat Chimeric.out.junction | awk '($1 ~/AB513134/ || $4 ~/AB513134/) && $1!="chrM" && $4!="chrM" && $7>0 && $8+$9<=5 {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | sort | uniq -c | sort -k1,1rn > Chimeric.out.junction.Count
    echo
    echo "Filter reads which map to HTLV-1 only"
        awk -F"\t" '($3=="AB513134"|| $7=="AB513134" || $1~/^@/) {print}' accepted_hits_sort_uniq.sam > htlv-1_accepted_hits_sort_uniq.sam
    echo
    echo "Read 1 - sense selection"
        awk '($2 ~/163/) || ($2 ~/83/) || ($1 ~/^@/) {print}' htlv-1_accepted_hits_sort_uniq.sam > output_read1_sense.sam
    echo
    echo "Read 1 - antisense selection"
        awk '($2 ~/99/) || ($2 ~/147/) || ($1 ~/^@/) {print}' htlv-1_accepted_hits_sort_uniq.sam > output_read1_antisense.sam
    echo
	echo "Analysis complete - displaying read counts"
	echo
	echo "Mapped reads after rmdup"
		samtools view -F 0x4 accepted_hits_sort_uniq.bam | cut -f 1 | sort | uniq | wc -l
	echo
	echo "Reads mapped to HTLV-1"
	echo
	    samtools view -Sb htlv-1_accepted_hits_sort_uniq.sam > htlv-1_accepted_hits_sort_uniq.bam
	echo
		samtools view -F 0x4 htlv-1_accepted_hits_sort_uniq.bam | cut -f 1 | sort | uniq | wc -l
 	echo
    echo "Read 1 - sense"
    echo
        samtools view -Sb output_read1_sense.sam > output_read1_sense.bam
    echo
        samtools view -F 0x4 output_read1_sense.bam | cut -f 1 | sort | uniq | wc -l
    echo
    echo "Read 1 - antisense"
        samtools view -Sb output_read1_antisense.sam > output_read1_antisense.bam
    echo
        samtools view -F 0x4 output_read1_antisense.bam | cut -f 1 | sort | uniq | wc -l
    echo
        echo
        echo "Congratulations! Data processing complete!"
        echo "Task completed on"
                date
        end=$(date +%s)
        duration=$(($end-$begin))
        echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
        echo
#
# ************ Main Script Ends Here ************ #
