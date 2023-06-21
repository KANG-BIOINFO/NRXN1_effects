library('biomaRt')
library('DESeq2')
library(readxl)
library(openxlsx)

setwd('/Volumes/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/rsem_transcripts/differential_v2_DEG/')

# prepare gene / transcript metadata for conversion
ensembl = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl', mirror = 'useast')
tps_meta = getBM(attributes = c(
  'ensembl_gene_id',
  'hgnc_symbol'), mart = ensembl)
tps_meta = tps_meta[!duplicated(tps_meta$ensembl_gene_id),]
rownames(tps_meta) = tps_meta$ensembl_gene_id

# convert for DE_isoforms (donor)
donor_res = readRDS('./DESeq2_result_donor.rds')
genes_sim = sapply(strsplit(rownames(donor_res), '\\.'), `[`, 1)
donor_res['gene_symbol'] = tps_meta[genes_sim, 'hgnc_symbol']
write.xlsx(as.data.frame(donor_res), 'DEGs_donor_combined.xlsx')

# convert for DE_isoforms (engineered)
donor_res = readRDS('./DESeq2_result_engineered.rds')
genes_sim = sapply(strsplit(rownames(donor_res), '\\.'), `[`, 1)
donor_res['gene_symbol'] = tps_meta[genes_sim, 'hgnc_symbol']
write.xlsx(as.data.frame(donor_res), 'DEGs_engineered.xlsx')

