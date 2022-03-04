#!/bin/bash

####################################################
#                                                  # 
#          kallisto for Paired-End reads           #
#				                   #  
#                                                  #
####################################################
#
# 
# ************ Main Script Begins Here ************ #
#
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
	    kallisto quant -i /reference_genome/kallisto/hg38_HTLV -l 75 -s 1 -t 10 -b 100 -o transcripts ../1_qc/clean_R_1.fastq ../1_qc/clean_R_2.fastq
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
