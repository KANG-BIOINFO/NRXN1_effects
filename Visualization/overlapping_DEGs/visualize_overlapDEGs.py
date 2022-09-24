import pandas as pd 
import numpy as np 
from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt 

df_overlap = pd.read_csv('num_pairwise_intersection_genes_donor.txt', sep = '\t', header = 0, index_col = 0)
df_meta = pd.read_csv('num_pairwise_intersection_genes_donor_metadata.txt', sep = '\t', header = 0, index_col = 0)

items = list(df_overlap.index.values)
n_items = len(items)
df_overlap = df_overlap.loc[:, items]
df_meta = df_meta.loc[items, :].copy()
df_meta['total_num_DEGs'] = np.max(df_overlap.values, axis = 1)

# hypergeometric tests
n_genes = 33538 # total number of genes
episilon = 1e-100
hyper_scores = pd.DataFrame(data = np.zeros((n_items, n_items)),
                            index = items, columns = items)
hyper_scores_log = pd.DataFrame(data = np.zeros((n_items, n_items)),
                                index = items, columns = items)                            
hyper_scores_corr = pd.DataFrame(data = np.zeros((n_items, n_items)),
                               index = items, columns = items)
hyper_scores_corr_log = pd.DataFrame(data = np.zeros((n_items, n_items)),
                                      index = items, columns = items)                                                           

df_p = pd.DataFrame(columns = ['i','j','pval'])
for i in range(n_items):
    for j in range(i, n_items):
        n_set1 = df_meta.loc[items[i], 'total_num_DEGs']
        n_set2 = df_meta.loc[items[j], 'total_num_DEGs']
        n_overlap = df_overlap.iloc[i, j]
        score = hypergeom(n_genes, n_set1, n_set2).sf(n_overlap - 1)
        hyper_scores.iloc[i, j] = hyper_scores.iloc[j, i] = score if score!=0 else episilon
        hyper_scores_log.iloc[i, j] = hyper_scores_log.iloc[j, i] = -np.log10(score)
        df_p.loc[df_p.shape[0], :] = [i, j, score]

# adjust p values
_, df_p['fdr'] = fdrcorrection(df_p['pval'])
for i in range(df_p.shape[0]):
    x, y, _, padj = df_p.iloc[i, :]
    # hyper_scores_corr.iloc[x, y] = padj
    # hyper_scores_corr_log.iloc[x, y] = -np.log10(padj)
    hyper_scores_corr.iloc[x, y] = hyper_scores_corr.iloc[y, x] = padj
    hyper_scores_corr_log.iloc[x, y] = hyper_scores_corr_log.iloc[y, x] = -np.log10(padj)

# format
hyper_scores_log.values[np.isinf(hyper_scores_log)] = 100
hyper_scores_log = hyper_scores_log.clip(0,100)

hyper_scores_corr_log.values[np.isinf(hyper_scores_corr_log)] = 100
hyper_scores_corr_log = hyper_scores_corr_log.clip(0,100)

df_meta.to_csv('num_pairwise_intersection_genes_donor_metadata_v2.txt',sep='\t')
hyper_scores.to_csv('hypergeometric_pvals_pairwise_intersection_genes_donor.txt', sep = '\t')
hyper_scores_log.to_csv('hypergeometric_logpvals_pairwise_intersection_genes_donor.txt', sep = '\t')
# hyper_scores_corr.to_csv('hypergeometric_fdr_pairwise_intersection_genes_donor.txt', sep = '\t')
# hyper_scores_corr_log.to_csv('hypergeometric_logfdr_pairwise_intersection_genes_donor.txt', sep = '\t')

hyper_scores_corr.to_csv('hypergeometric_fdr_pairwise_intersection_genes_donor_v2.txt', sep = '\t') # make symmetric instead of triangle shape
hyper_scores_corr_log.to_csv('hypergeometric_logfdr_pairwise_intersection_genes_donor_v2.txt', sep = '\t')
