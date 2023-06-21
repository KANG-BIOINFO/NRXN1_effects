library(tximport)
library(DESeq2)
library(stringr)
library(openxlsx)

setwd('/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/rsem_transcripts/differential_v2')
source_folder = '/data/aronow/Kang/induced_neuron/ChangHui_RNAseq_fastq/isoform_RSEM/rsem_transcripts/transcripts/'

# load data
df_samples = read.table('../samples_temporary.txt', sep = '\t', header = 1, row.names = 1)
df_donors = df_samples[df_samples$iN.subtype=='direct',]
df_eng = df_samples[df_samples$iN.subtype=='engineered',]

# donor
samples = rownames(df_donors)
isoform_files = paste0('../transcripts/', samples, '.isoforms.results')
names(isoform_files) = samples
txi <- tximport(files = isoform_files, type = "rsem", txIn = T, txOut = T)
df_donors$NRXN1.genome = factor(df_donors$NRXN1.genome, levels = c('WT', 'HET'))
dds <- DESeqDataSetFromTximport(txi, colData = df_donors, design = ~CellLine + NRXN1.genome)
dds <- DESeq(dds)
res = results(dds)
saveRDS(dds, 'DESeq2_object_donor.rds')
saveRDS(res, 'DESeq2_result_donor.rds')
write.xlsx(res, 'DE_isoforms_donor.xlsx')

# eng
samples = rownames(df_eng)
isoform_files = paste0('../transcripts/', samples, '.isoforms.results')
names(isoform_files) = samples
txi <- tximport(files = isoform_files, type = "rsem", txIn = T, txOut = T)
df_eng$NRXN1.genome = factor(df_eng$NRXN1.genome, levels = c('WT', 'HET'))
dds <- DESeqDataSetFromTximport(txi, colData = df_eng, design = ~ NRXN1.genome)
dds <- DESeq(dds)
res = results(dds)
saveRDS(dds, 'DESeq2_object_engineered.rds')
saveRDS(res, 'DESeq2_result_engineered.rds')
write.xlsx(res, 'DE_isoforms_engineered.xlsx')
