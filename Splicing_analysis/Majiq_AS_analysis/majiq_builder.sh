#!/bin/bash

#BSUB -M 128000
#BSUB -W 12:00
#BSUB -n 12
#BSUB -R "rusage[mem=128000] span[hosts=1]"
#BSUB -J majiq_builder
#BSUB -o logging/majiq_builder.out
#BSUB -e logging/majiq_builder.err

module load anaconda3
source activate /data/aronow/Kang/py_envir/splice2/

# build
n_jobs=8
gff_file='/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/run_majiq/source_data/gff_file/gencode.v42.annotation.gff3'
output_folder='/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/run_majiq/output/spliceGraph/'
majiq build  -c settings_file.ini $gff_file -o $output_folder  -j $n_jobs