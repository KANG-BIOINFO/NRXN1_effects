import os 
import math
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib as mpl
import seaborn

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams.update({'font.size': 24})

cmap = mpl.cm.get_cmap('tab20')

# df_genes = pd.read_excel('SCZ_schema-PGC3_gene_list_v2.xlsx')
df_selected = pd.read_excel('selection_genes.xlsx')
df_selected.drop_duplicates(subset = ['gene'], inplace = True)
df_example = pd.read_csv('/Users/kang/Dropbox/scRNA-seq/NRXN1 single cell RNA-seq with Bruce (1)/single cell analysis/DEG/Jul4_time_samplegroup_celltype_Case_vs_Ctrl_4_timpoints/d22_Cycling Dorsal NEC_Donor-NRXN1del_vs_Donor-Ctrl.txt', sep = '\t', header = 0, index_col=0)

overlap_genes = np.intersect1d(df_selected['gene'], df_example.index.values)
df_selected.set_index('gene', inplace = True)
available_genes = [i for i in df_selected.index.values if i in overlap_genes]
df_selected = df_selected.loc[available_genes, :].copy()

# source_order = [('Donor-NRXN1del_vs_Donor-Ctrl', 'd22'),
#                 ('Donor-NRXN1del_vs_Donor-Ctrl', 'd50'),
#                 ('Donor-NRXN1del_vs_Donor-Ctrl', 'd101'),
#                 ('Engineered-Cre_vs_Engineered-Flp', 'd23'),
#                 ('Engineered-Cre_vs_Engineered-Flp', 'd50'),
#                 ('Engineered-Cre_vs_Engineered-Flp', 'd112'),]

# only draw plots for d101 donor and d112 eng
source_order = [('Donor-NRXN1del_vs_Donor-Ctrl', 'd101'),
                ('Engineered-Cre_vs_Engineered-Flp', 'd112')]
cellclass_ordered_v2 = ["Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",\
                         "Astroglia", "Glia cell", "oRG", "Intermediate cells", "Choroid Plexus",\
                         "OPC", "IPC1","IPC2","IPC3",\
                         "IN1","IN2","IN3","IN4","IN5","IN6","IN7",\
                         "CN1","CN2","CN3","CN4","CN5"]
folder_name = '/Users/kang/Dropbox/scRNA-seq/NRXN1 single cell RNA-seq with Bruce (1)/single cell analysis/DEG/Oct12_time_samplegroup_celltype_Case_vs_Ctrl/'

df_results = df_selected.copy()
fig, axes = plt.subplots(figsize = (25, 30), nrows = 1, ncols = 2, sharey = True)
for idx, source in enumerate(source_order):
    organoid_type, timepoint = source[0], source[1]
    i = 0
    xlabels = []
    x_all = []
    y_all = []
    s_all = []
    c_all = []

    for cellclass in cellclass_ordered_v2:
        file_name = folder_name + timepoint + '_' + cellclass + '_' + organoid_type + '.txt'
        if os.path.exists(file_name):
            df_deg = pd.read_csv(file_name, sep = '\t', header = 0, index_col = 0)
            df_deg = df_deg.loc[available_genes, :].copy()
        
            scores = [0 if j > 0.05 else -math.log10(j) for j in df_deg['pvals_adj']]
            logfcs = list(df_deg['logfoldchanges'])

            # axes[idx].scatter(
            #     x = [i] * len(available_genes), 
            #     y = np.arange(len(available_genes)), 
            #     s = np.array(scores) * 25, 
            #     c = [cmap(i / 26)] * len(available_genes),
            #     label = cellclass
            # )
            for k, score, logfc in zip(np.arange(len(available_genes)), scores, logfcs):
                if score != 0: 
                    axes[idx].scatter(x = i, y = k, s = score * 25, c = cmap(i / 26))#, c = logfcs)
                    x_all.append(i); y_all.append(k), s_all.append(score); c_all.append(cmap(i / 26))
            i += 1
            xlabels.append(cellclass)

    legend_sizes = np.sort(s_all)[::len(s_all) // 4][-3:]
    indices = [np.where(s_all == v)[0][0] for v in legend_sizes]

    for l in indices:
        axes[idx].scatter(x = x_all[l], y = y_all[l], s = 25 * s_all[l], c = c_all[l], label = '{:.1f}'.format(s_all[l]))

    axes[idx].set_xticks(np.arange(len(xlabels)))
    axes[idx].set_xticklabels(xlabels, rotation = 90)
    axes[idx].set_xlim(-0.5, len(xlabels))
    axes[idx].set_xlabel(organoid_type + '-' + timepoint)
    axes[idx].legend(scatterpoints = 1, title = 'DE Significance\n (-log10(pval_adj))')
    axes[idx].grid(linewidth = 0.5, alpha = 0.5)
    
# plt.xticks(rotation = 90)
plt.yticks(np.arange(len(available_genes)), available_genes)
plt.ylim(-0.5, len(available_genes) + 0.5)
plt.tight_layout()

fig.savefig('selected_genes/dotplot_d101_donor_d112_eng_v2.pdf', bbox_inches = 'tight')

