#!/bin/bash

#BSUB -M 128000
#BSUB -W 12:00
#BSUB -n 6
#BSUB -R "rusage[mem=128000] span[hosts=1]"
#BSUB -J rsem_de_isoforms
#BSUB -o rsem_de_isoforms.out
#BSUB -e rsem_de_isoforms.err

source_folder='transcripts/'
output_folder='differential/'
module load RSEM/1.3.0

# cell line A
rsem-generate-data-matrix ${source_folder}'19N-55429417.isoforms.results' \
                            ${source_folder}'26N-55424500.isoforms.results' ${source_folder}'N9540_19-55430413.isoforms.results' \
                            ${source_folder}'19C-55431414.isoforms.results' ${source_folder}'26C-55439384.isoforms.results' \
                            ${source_folder}'C2463_19-55427424.isoforms.results' > ${output_folder}'IsoMat_LineA.txt'

# rsem-run-ebseq ${output_folder}'GeneMat_LineA.txt' 3,3 ${output_folder}'GeneMat_LineA.results'
# rsem-control-fdr ${output_folder}'GeneMat_LineA.results' 0.05 ${output_folder}'GeneMat_LineA.de.txt'

# cell line B
rsem-generate-data-matrix ${source_folder}'N1-Old-TAGCTT_S10.isoforms.results' \
                            ${source_folder}'N2-Old-GGCTAC_S11.isoforms.results' ${source_folder}'N3-Old-CTTGTA_S12.isoforms.results' \
                            ${source_folder}'C1-Old-CAGATC_S7.isoforms.results' ${source_folder}'C2-Old-ACTTGA_S8.isoforms.results' \
                            ${source_folder}'C3-Old-GATCAG_S9.isoforms.results' > ${output_folder}'IsoMat_LineB.txt'

# rsem-run-ebseq ${output_folder}'GeneMat_LineB.txt' 3,3 ${output_folder}'GeneMat_LineB.results'
# rsem-control-fdr ${output_folder}'GeneMat_LineB.results' 0.05 ${output_folder}'GeneMat_LineB.de.txt'

# cell line C
rsem-generate-data-matrix ${source_folder}'N3320_23-55429418.isoforms.results' \
                            ${source_folder}'N3320_4-55437387.isoforms.results' ${source_folder}'N3320_7-55421553.isoforms.results' \
                            ${source_folder}'C3141_23-55438387.isoforms.results' ${source_folder}'C3141_4-55428435.isoforms.results' \
                            ${source_folder}'C3141_7-55431415.isoforms.results' > ${output_folder}'IsoMat_LineC.txt'

# rsem-run-ebseq ${output_folder}'GeneMat_LineC.txt' 3,3 ${output_folder}'GeneMat_LineC.results'
# rsem-control-fdr ${output_folder}'GeneMat_LineC.results' 0.05 ${output_folder}'GeneMat_LineC.de.txt'

# cell line D
rsem-generate-data-matrix ${source_folder}'1C-New-TGACCA_S4.isoforms.results' \
                            ${source_folder}'2C-New-ACAGTG_S5.isoforms.results' ${source_folder}'3C-New-GCCAAT_S6.isoforms.results' \
                            ${source_folder}'1F-New-ATCACG_S1.isoforms.results' ${source_folder}'2F-New-CGATGT_S2.isoforms.results' \
                            ${source_folder}'3F-New-TTAGGC_S3.isoforms.results' > ${output_folder}'IsoMat_LineD.txt'
