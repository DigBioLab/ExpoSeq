#!/bin/bash

fastq_directory=$1 	
path_to_mixcr=$2
save_dir=$PWD
paired_end_sequencing=False		#
threads=1	#
method="milab-human-tcr-dna-multiplex-cdr3"
trim_div_by=3 #
trim_min_count=3 #

python "mixcr_cl.py" $fastq_directory $path_to_mixcr $save_dir $paired_end_sequencing $threads $method $trim_div_by $trim_min_count