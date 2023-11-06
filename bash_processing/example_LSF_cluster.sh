#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J darian2
### -- ask for number of cores (default: 1) -- 
#BSUB -n 5
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=10GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
#BSUB -M 10GB
### -- set walltime limit: hh:mm -- 
#BSUB -W 2:00 
### -- set the email address --
#BSUB -u s220672@dtu.dk
### -- send notification at start -- 
#BSUB -B 
### -- send notification at completion -- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o job_output/darian2.out 
#BSUB -e job_output/darian2.err 


module load python3/3.9.11
module swap openblas/0.3.7
module swap ffmpeg/4.2.2   
module load numpy/1.22.3-python-3.9.11-openblas-0.3.19
module load pandas/1.0.3-python-3.8.2
module load openblas/0.3.19

script_name="~/ExpoSeq/bash_processing/mixcr_cl.py" 
mixcr_path="~/mixcr/mixcr.jar"
sample_name="test_darian_exome_full_length"
fastq_forward="~/fastq_files/darian"
#fastq_reverse="~/fastq_files/FastqFiles_reverse_test"
threads=5
java_heap_size=3000
add_arguments="--species hsa"
method="ampliseq-tcrb-plus-full-length"


python3 $script_name $mixcr_path $sample_name $fastq_forward --threads $threads --java_heap_size $java_heap_size --method $method --add_arguments "--species hsa" --mixcr_chain "IG"