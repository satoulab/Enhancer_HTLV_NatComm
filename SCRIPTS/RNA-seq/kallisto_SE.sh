#!/bin/bash

####################################################
#                                                  # 
#          kallisto for Single-end reads	   #
#  				                   #  
#                                                  #
####################################################

# ************ Main Script Begins Here ************ #

	echo
	echo "Begin analysis"
	echo "Task started at"
                date
        begin=$(date +%s)
    	echo
# Make transcript quantification directory
        echo "Create directory 3_transcript_quant & begin quantification in new directory"
		mkdir 3_transcript_quant
                cd 3_transcript_quant

# Transcript quantification
	echo "run kallisto"
	echo
	    	kallisto quant --single -i /reference_genome/kallisto/HTLV/hg38_HTLV_cDNA -l 75 -s 1 -t 10 -b 100 -o transcripts ../1_qc/R1_clean.fastq
		cd ../1_qc
		gzip R1_clean.fastq
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
