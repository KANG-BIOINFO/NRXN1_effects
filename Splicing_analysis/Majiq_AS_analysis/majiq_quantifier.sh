#!/bin/bash

#BSUB -M 128000
#BSUB -W 12:00
#BSUB -n 12
#BSUB -R "rusage[mem=128000] span[hosts=1]"
#BSUB -J majiq_quantifier
#BSUB -o logging/majiq_quantifier.out
#BSUB -e logging/majiq_quantifier.err

module load anaconda3
source activate /data/aronow/Kang/py_envir/splice/

# build
input_folder='/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/run_majiq/output/spliceGraph/'
output_folder='/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/run_majiq/output/heterogen/'
majiq heterogen -o output/heterogen \
    -n A_HET A_WT \
    -grp1 ${input_folder}{19N-55429417.sorted,26N-55424500.sorted,N9540_19-55430413.sorted}.majiq \
    -grp2 ${input_folder}{19C-55431414.sorted,26C-55439384.sorted,C2463_19-55427424.sorted}.majiq \

majiq heterogen -o output/heterogen \
    -n B_HET B_WT \
    -grp1 ${input_folder}{N1-Old-TAGCTT_S10.sorted,N2-Old-GGCTAC_S11.sorted,N3-Old-CTTGTA_S12.sorted}.majiq \
    -grp2 ${input_folder}{C1-Old-CAGATC_S7.sorted,C2-Old-ACTTGA_S8.sorted,C3-Old-GATCAG_S9.sorted}.majiq \

majiq heterogen -o output/heterogen \
    -n C_HET C_WT \
    -grp1 ${input_folder}{N3320_23-55429418.sorted,N3320_4-55437387.sorted,N3320_7-55421553.sorted}.majiq \
    -grp2 ${input_folder}{C3141_23-55438387.sorted,C3141_4-55428435.sorted,C3141_7-55431415.sorted}.majiq \

majiq heterogen -o output/heterogen \
    -n D_HET D_WT \
    -grp1 ${input_folder}{1C-New-TGACCA_S4.sorted,2C-New-ACAGTG_S5.sorted,3C-New-GCCAAT_S6.sorted}.majiq \
    -grp2 ${input_folder}{1F-New-ATCACG_S1.sorted,2F-New-CGATGT_S2.sorted,3F-New-TTAGGC_S3.sorted}.majiq \