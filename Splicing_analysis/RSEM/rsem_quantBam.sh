#!/bin/bash

#BSUB -M 500000
#BSUB -W 12:00
#BSUB -n 12
#BSUB -R "rusage[mem=500000] span[hosts=1]"
#BSUB -J rsem_quant
#BSUB -o rsem_quant.out
#BSUB -e rsem_quant.err

module load RSEM/1.3.0

BAM_folder='/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/star_align/BAM_transcripts/'
BAM_file='Aligned.toTranscriptome.out.bam'
output_folder='/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/rsem_transcripts/transcripts/'
reference='/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/rsem_transcripts/rsem_reference/RSEM_STAR_Ref/hg38'
cd $BAM_folder
samples=$(find . -maxdepth 1 -name "*-*")

foo() {
    local sample_dummy=$1
    sample=${sample_dummy:2:${#sample_dummy}}
    cd $BAM_folder$sample'/'
    
    rsem-calculate-expression -p 12 \
    --bam \
    --paired-end \
    --estimate-rspd \
    --no-bam-output \
    $BAM_file \
    $reference $output_folder$sample
}

for sample_dummy in $samples; do
    foo "$sample_dummy" &
done
