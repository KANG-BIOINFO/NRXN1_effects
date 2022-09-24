import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# settings 
cellclass_ordered_v2 = ["Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                         "Astroglia", "Glia cell", "oRG", "Choroid Plexus", "Intermediate cells",
                         "OPC", "IPC1","IPC2","IPC3",
                         "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                         "CN1","CN2","CN3","CN4","CN5"]

# load data
DEG_original_folder = '/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/DEGs/survey/DEG_num/'
def draw_heatmap_up_n_dn(deg_files, labels, day = ['d112']):
    # draw up / down / total for d112 in the same heatmap
    df_combined = pd.DataFrame()
    for deg_file, label in zip(deg_files, labels):
        df = pd.read_csv(DEG_original_folder + deg_file, sep = '\t', 
                        header = 0, index_col = 0, keep_default_na = True)
        df = df.fillna(0)
        df = df.loc[cellclass_ordered_v2, : ].copy()
        df = df.loc[:, day].copy()
        df.columns = [label]
        df_combined = pd.concat([df_combined, df], axis = 1)

    fig, ax = plt.subplots(figsize = (3.5, 4.8))
    # sns.heatmap(df, cmap = 'YlGnBu', ax = ax)
    sns.heatmap(df_combined, cmap = 'Greens', ax = ax)
    fig.savefig('figures/DEG_num_Eng_d112_all.pdf', bbox_inches = 'tight')

deg_files = [
    'DEG_num_down_engineered.txt',
    'DEG_num_up_engineered.txt',
    'DEG_num_total_engineered.txt'
]
labels = ['Down', 'Up', 'Total']

# draw_heatmap_up_n_dn(deg_files, labels, day = ['d112'])

# draw DEG numbers for filtered table
df_deg = pd.read_excel('/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/DEGs/survey/DEGs_source_time_celltype_ribosomalTotallyRemoved.xlsx')
df_deg = df_deg.loc[(df_deg['source'] == 'Engineered-Cre') & (df_deg['timepoint'] == 'd112'), : ].copy()

df_deg_num = pd.DataFrame(
    data = np.zeros((len(cellclass_ordered_v2), 3)),
    index = cellclass_ordered_v2,
    columns = ['Down', 'Up', 'Total']
)

for cellclass in cellclass_ordered_v2:
    df_deg_sub = df_deg.loc[df_deg['cell_type'] == cellclass, : ].copy()
    df_deg_num.loc[cellclass, 'Down'] = df_deg_sub.loc[df_deg_sub['logfoldchanges'] > 0, :].shape[0]
    df_deg_num.loc[cellclass, 'Up'] = df_deg_sub.loc[df_deg_sub['logfoldchanges'] < 0, :].shape[0]
    df_deg_num.loc[cellclass, 'Total'] = df_deg_sub.shape[0]

def draw_heatmap_v2(df):
    df = df.fillna(0)
    df = df.loc[cellclass_ordered_v2, : ].copy()

    fig, ax = plt.subplots(figsize = (3.5, 4.8))
    # sns.heatmap(df, cmap = 'YlGnBu', ax = ax)
    sns.heatmap(df, cmap = 'Greens', ax = ax)
    fig.savefig('figures/DEG_num_Eng_d112_neurogenesis_genes.pdf', bbox_inches = 'tight')

draw_heatmap_v2(df_deg_num)
