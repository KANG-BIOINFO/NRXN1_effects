#!/bin/bash

#BSUB -M 64000
#BSUB -W 12:00
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -J STAR_align_batch1
#BSUB -o STAR_align_batch1.out
#BSUB -e STAR_align_batch1.err

module load STAR/2.7.9

genomeDir_='/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/star_reference/GRCh38_p13_index'
num_threads=8
parent_folder='/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/Bulk_RNAseq/fastq_batch1/'

cd $parent_folder
samples=$(find . -maxdepth 1 -name "*-*")
for sample_dummy in $samples; do
    sample=${sample_dummy:2:${#sample_dummy}}
    mkdir '/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/star_align/BAM/'$sample
    outputFile='/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/star_align/BAM/'$sample'/'

    sample_folder=${parent_folder}$sample'/'
    cd $sample_folder
    R1_files=$(find . -maxdepth 1 -name "*R1*")
    R2_files=$(find . -maxdepth 1 -name "*R2*")
    R1_files=$(sort <<<"${R1_files[*]}")
    R2_files=$(sort <<<"${R2_files[*]}")

    R1_files_concat=()
    for i in $R1_files; do
        R1_files_concat+=$i','
    done
    R1_files_concat=${R1_files_concat::-1}

    R2_files_concat=()
    for i in $R2_files; do
        R2_files_concat+=$i','
    done
    R2_files_concat=${R2_files_concat::-1}

    input_files=$R1_files_concat' '$R2_files_concat
    echo $input_files

    STAR --genomeDir $genomeDir_\
    --runThreadN $num_threads \
    --readFilesCommand zcat \
    --readFilesIn $input_files \
    --outFileNamePrefix $outputFile \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard 
done
