# Setup
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz

fdr_vals = [0.05, 0.10, 0.20, 0.30, 0.40]
cellclass_ordered_v2 =  ["Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                         "Astroglia", "Glia cell", "oRG", "Intermediate cells", "Choroid Plexus",
                         "OPC", "IPC1","IPC2","IPC3",
                         "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                         "CN1","CN2","CN3","CN4","CN5"]

cellclass_ordered =  ["Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                      "Astroglia", "Glia cell", "oRG", "Intermediate cells", "Choroid Plexus",
                      "Mesenchymal cell", "Microglia",
                      "OPC", "IPC1","IPC2","IPC3",
                      "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                      "CN1","CN2","CN3","CN4","CN5"]                         

def run_sccoda_(df_donor, timepoint, cellclass_list, ref_type, suffix = ''):
    # prepare input table
    df_donor_subset = df_donor.loc[df_donor['cellclass_draft_v2'].isin(cellclass_list), : ].copy()
    df_donor_subset = df_donor_subset.loc[df_donor_subset['timepoint_v2'] == timepoint,:].copy()
    df_counts = df_donor_subset.groupby(['sample_v2','cellclass_draft_v2']).size()
    df_counts = df_counts.unstack()
    df_counts.values[np.isnan(df_counts.values)] = 0

    df_counts.to_csv('test_cell_counts.txt', sep = '\t')
    df_counts = pd.read_csv('test_cell_counts.txt', sep = '\t', header = 0, keep_default_na = True)

    data_all = dat.from_pandas(df_counts, covariate_columns=["sample_v2"])
    data_all.obs['condition'] = ['Control' if i.startswith('Control') else 'NRXN1_del' for i in data_all.obs['sample_v2']]

    # visualize
    viz.boxplots(data_all, feature_name = 'condition', figsize = (12, 8))
    plt.savefig('figures/boxplot_' + timepoint + suffix + '.pdf', bbox_inches = 'tight')

    # model setup
    model = mod.CompositionalAnalysis(data_all, formula = 'condition', reference_cell_type = ref_type)
    sim_results = model.sample_hmc()

    for fdr_val in fdr_vals:
        sim_results.set_fdr(est_fdr = fdr_val)
        sim_results.effect_df.to_csv('results/effects_' + timepoint + '_fdr_' + str(fdr_val) + suffix + '.txt', sep = '\t')

    sim_results.save('results/result_' + timepoint + suffix)

if __name__ == '__main__':
    # load data
    df_all = pd.read_csv('/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/ctype_v6.txt', sep = '\t', header = 0, index_col = 0)
    df_eng = df_all.loc[df_all['source'] == 'Engineered', : ].copy()
    df_eng['sample_v2'] = [(i + '_' + j) for (i, j) in zip(df_eng['genotype'], df_eng['sample'])]

    timepoints = ['d023', 'd050', 'd112']
    ref_types = ['Cycling Dorsal NEC', 'Cycling Dorsal NEC', 'Cycling Dorsal NEC']
    for i in range(3):
        timepoint = timepoints[i]
        ref_type = ref_types[i]
        run_sccoda_(df_eng, timepoint, cellclass_ordered, ref_type, suffix = '_Engineered_cellclass_1')
        run_sccoda_(df_eng, timepoint, cellclass_ordered_v2, ref_type, suffix = '_Engineered_cellclass_2')
    