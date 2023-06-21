#!/bin/bash

#BSUB -M 128000
#BSUB -W 12:00
#BSUB -n 6
#BSUB -R "rusage[mem=128000] span[hosts=1]"
#BSUB -J rsem_de_transcripts
#BSUB -o rsem_de_transcripts.out
#BSUB -e rsem_de_transcripts.err

source_folder='transcripts/'
output_folder='differential/'
module load RSEM/1.3.0

# cell line A
rsem-generate-data-matrix ${source_folder}'19N-55429417.genes.results' \
                            ${source_folder}'26N-55424500.genes.results' ${source_folder}'N9540_19-55430413.genes.results' \
                            ${source_folder}'19C-55431414.genes.results' ${source_folder}'26C-55439384.genes.results' \
                            ${source_folder}'C2463_19-55427424.genes.results' > ${output_folder}'GeneMat_LineA.txt'

# rsem-run-ebseq ${output_folder}'GeneMat_LineA.txt' 3,3 ${output_folder}'GeneMat_LineA.results'
# rsem-control-fdr ${output_folder}'GeneMat_LineA.results' 0.05 ${output_folder}'GeneMat_LineA.de.txt'

# cell line B
rsem-generate-data-matrix ${source_folder}'N1-Old-TAGCTT_S10.genes.results' \
                            ${source_folder}'N2-Old-GGCTAC_S11.genes.results' ${source_folder}'N3-Old-CTTGTA_S12.genes.results' \
                            ${source_folder}'C1-Old-CAGATC_S7.genes.results' ${source_folder}'C2-Old-ACTTGA_S8.genes.results' \
                            ${source_folder}'C3-Old-GATCAG_S9.genes.results' > ${output_folder}'GeneMat_LineB.txt'

# rsem-run-ebseq ${output_folder}'GeneMat_LineB.txt' 3,3 ${output_folder}'GeneMat_LineB.results'
# rsem-control-fdr ${output_folder}'GeneMat_LineB.results' 0.05 ${output_folder}'GeneMat_LineB.de.txt'

# cell line C
rsem-generate-data-matrix ${source_folder}'N3320_23-55429418.genes.results' \
                            ${source_folder}'N3320_4-55437387.genes.results' ${source_folder}'N3320_7-55421553.genes.results' \
                            ${source_folder}'C3141_23-55438387.genes.results' ${source_folder}'C3141_4-55428435.genes.results' \
                            ${source_folder}'C3141_7-55431415.genes.results' > ${output_folder}'GeneMat_LineC.txt'

# rsem-run-ebseq ${output_folder}'GeneMat_LineC.txt' 3,3 ${output_folder}'GeneMat_LineC.results'
# rsem-control-fdr ${output_folder}'GeneMat_LineC.results' 0.05 ${output_folder}'GeneMat_LineC.de.txt'

# cell line D
rsem-generate-data-matrix ${source_folder}'1C-New-TGACCA_S4.genes.results' \
                            ${source_folder}'2C-New-ACAGTG_S5.genes.results' ${source_folder}'3C-New-GCCAAT_S6.genes.results' \
                            ${source_folder}'1F-New-ATCACG_S1.genes.results' ${source_folder}'2F-New-CGATGT_S2.genes.results' \
                            ${source_folder}'3F-New-TTAGGC_S3.genes.results' > ${output_folder}'GeneMat_LineD.txt'

# rsem-run-ebseq ${output_folder}'GeneMat_LineD.txt' 3,3 ${output_folder}'GeneMat_LineD.results'
# rsem-control-fdr ${output_folder}'GeneMat_LineD.results' 0.05 ${output_folder}'GeneMat_LineD.de.txt'
