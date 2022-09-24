import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt 
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

import seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 12

# load and extract info
df_enrich = pd.read_excel('disease_enrichment_selection_reformat.xlsx')

index_names = list(df_enrich['Unnamed: 0'])
index_cats = [i.split('_')[0] for i in index_names]
index_des = ['_'.join(i.split('_')[1:]) for i in index_names]

df_enrich.set_index('Unnamed: 0', inplace = True)
df_enrich.index.name = 'enrichment_terms'

col_names = df_enrich.columns.tolist()
col_time = [i.split('_')[0] for i in col_names]
col_cells = [i.split('_')[1] for i in col_names]
col_organoids = ['Donor' if i == 'd101' else 'Engineered' for i in col_names]
col_regulations = [i.split('_')[-1].split('-regulation')[0] for i in col_names]

df_enrich.columns = col_cells

# draw plot
thred_1 = -np.log10(0.05)
thred_2 = -np.log10(0.01)
thred_3 = -np.log10(0.001)

vlag = sns.color_palette("vlag", as_cmap=True)
vlag_blue = vlag(np.linspace(0, 1, 102))[:51]
vlag_red = vlag(np.linspace(0, 1, 410))[205:]
newcolors = np.vstack((vlag_blue, vlag_red))
newcmp = ListedColormap(newcolors, name = 'redblue')

df_enrich = df_enrich.T.copy()

fig, ax = plt.subplots(figsize = (20, 4))
sns.heatmap(
    df_enrich, 
    ax = ax, 
    linewidths = 0.3,
    linecolor = 'black',
    vmin = 0,
    vmax = 10,
    cmap = newcmp
)

for i in range(df_enrich.shape[0]):
    for j in range(df_enrich.shape[1]):
        val = df_enrich.iloc[i, j]
        symbol = ''
        if val > thred_1 and val <= thred_2:
            symbol = '*' 
            ax.text(j + 0.40, i + 1 - 0.05, s = symbol)
        elif val > thred_2 and val <= thred_3:
            symbol = '**' 
            ax.text(j + 0.35, i + 1 - 0.05, s = symbol)
        elif val > thred_3:
            symbol = '***' 
            ax.text(j + 0.25, i + 1 - 0.25, s = symbol)
        

ax.set_xlabel('')
ax.xaxis.set_ticks_position('top') 
plt.xticks(rotation = 90, ha = 'right')
ax.set_ylabel('')
ax.yaxis.set_ticks_position('right')
plt.yticks(rotation = 360, va = 'top')
fig.savefig('enrichment_heatmap.pdf', bbox_inches = 'tight')
