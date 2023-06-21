library('biomaRt')
library('DESeq2')
library(readxl)
library(openxlsx)

setwd('/Volumes/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/rsem_transcripts/differential_v2/')

# prepare gene / transcript metadata for conversion
ensembl = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl', mirror = 'useast')
tps_meta = getBM(attributes = c(
  'ensembl_gene_id',
  'ensembl_transcript_id',
  'hgnc_symbol'), mart = ensembl)
tps_meta = tps_meta[!duplicated(tps_meta$ensembl_transcript_id),]
rownames(tps_meta) = tps_meta$ensembl_transcript_id

# convert for DE_isoforms (donor)
donor_res = readRDS('./DESeq2_result_donor.rds')
transcripts_sim = sapply(strsplit(rownames(donor_res), '\\.'), `[`, 1)
donor_res['transcript_id'] = transcripts_sim
donor_res['gene_id'] = tps_meta[transcripts_sim, 'ensembl_gene_id']
donor_res['gene_symbol'] = tps_meta[transcripts_sim, 'hgnc_symbol']
write.xlsx(donor_res, 'DE_isoforms_donor.xlsx')

# convert for DE_isoforms (engineered)
donor_res = readRDS('./DESeq2_result_engineered.rds')
transcripts_sim = sapply(strsplit(rownames(donor_res), '\\.'), `[`, 1)
donor_res['transcript_id'] = transcripts_sim
donor_res['gene_id'] = tps_meta[transcripts_sim, 'ensembl_gene_id']
donor_res['gene_symbol'] = tps_meta[transcripts_sim, 'hgnc_symbol']
write.xlsx(donor_res, 'DE_isoforms_engineered.xlsx')

