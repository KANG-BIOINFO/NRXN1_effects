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
def draw_heatmap(dataframe_name, label):
    df = pd.read_csv(DEG_original_folder + dataframe_name, sep = '\t', 
                     header = 0, index_col = 0, keep_default_na = True)
    df = df.fillna(0)
    df = df.loc[cellclass_ordered_v2, : ].copy()

    fig, ax = plt.subplots(figsize = (3.5, 4.8))
    # sns.heatmap(df, cmap = 'YlGnBu', ax = ax)
    sns.heatmap(df, cmap = 'Greens', ax = ax)
    fig.savefig('figures/' + label +'.pdf', bbox_inches = 'tight')

deg_files = [
    'DEG_num_down_donor.txt',
    'DEG_num_up_donor.txt',
    'DEG_num_total_donor.txt'
]
labels = ['DEG_num_down_donor', 'DEG_num_up_donor', 'DEG_num_total_donor']

for i, j in zip(deg_files, labels):
    draw_heatmap(i, j)

# draw DEG numbers for filtered table
df_deg = pd.read_excel('/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/DEGs/survey/DEGs_source_time_celltype_ribosomalTotallyRemoved.xlsx')
df_deg = df_deg.loc[df_deg['source'] == 'Donor-NRXN1del', : ].copy()

timepoints = ['d22', 'd50', 'd101']
df_deg_up = pd.DataFrame(
    data = np.zeros((len(cellclass_ordered_v2), len(timepoints))),
    index = cellclass_ordered_v2,
    columns = timepoints
)
df_deg_dn = df_deg_up.copy(); df_deg_total = df_deg_up.copy()

for cellclass in cellclass_ordered_v2:
    for timepoint in timepoints:
        df_deg_sub = df_deg.loc[(df_deg['cell_type'] == cellclass) & (df_deg['timepoint'] == timepoint), : ].copy()
        df_deg_up.loc[cellclass, timepoint] = df_deg_sub.loc[df_deg_sub['logfoldchanges'] > 0, :].shape[0]
        df_deg_dn.loc[cellclass, timepoint] = df_deg_sub.loc[df_deg_sub['logfoldchanges'] < 0, :].shape[0]
        df_deg_total.loc[cellclass, timepoint] = df_deg_sub.shape[0]

def draw_heatmap_v2(df, label):
    df = df.fillna(0)
    df = df.loc[cellclass_ordered_v2, : ].copy()

    fig, ax = plt.subplots(figsize = (3.5, 4.8))
    # sns.heatmap(df, cmap = 'YlGnBu', ax = ax)
    sns.heatmap(df, cmap = 'Greens', ax = ax)
    fig.savefig('figures/' + label +'.pdf', bbox_inches = 'tight')

deg_dfs = [
    df_deg_dn,
    df_deg_up,
    df_deg_total
]
labels = ['DEG_num_down_donor_filtered', 'DEG_num_up_donor_filtered', 'DEG_num_total_donor_filtered']

for i, j in zip(deg_dfs, labels):
    draw_heatmap_v2(i, j)
