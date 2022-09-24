import scanpy as sc 
import numpy as np 
import pandas as pd 
import sys

def format_DEGs(adata):
    keys = ["names","scores","logfoldchanges","pvals","pvals_adj","pts"]
    for i,key in enumerate(keys):
        a = pd.DataFrame(adata.uns["rank_genes_groups"][key]) # transfer to data frame
        b = pd.DataFrame(a.values.T.reshape(1,a.shape[0]*a.shape[1]).T) # reformat the data frame to one column
           
        if i == 0:
            b.columns = [key] # rename the column name
            b["Status"] = sorted(list(a.columns)*a.shape[0]) # add Status annotation
            b.set_index([key],inplace=True)
            b_merged = b
        else:
            if key in ["pts"]:
                pts_all = []
                for cell_group in np.unique(b_merged["Status"]):
                    genes = b_merged.loc[b_merged["Status"] == cell_group,:].index.values
                    pts_all = pts_all + list(a.loc[genes, cell_group])
                b_merged[key] = pts_all
            else:
                b_merged[key] = list(b[0])
        
    return b_merged


if __name__ == "__main__":
    # Example: /data/aronow/Kang/py_envir/aipy/bin/python3.7 DE_markers.py "/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/combined_all_organoids_v3.h5ad" "integrated_snn_res.1" "13" "3"

    # 0. initialization
    adata = sc.read("/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/combined_all_organoids_v4.h5ad")
    method = "wilcoxon" # [‘logreg’, ‘t-test’, ‘wilcoxon’, ‘t-test_overestim_var’]

    meta = adata.obs.copy()
    n_genes = adata.shape[1]
    
    L_time = "timepoint"
    L_sampleG = "sample_group"
    L_cellT = "cellclass_draft_v2"

    group_engineered = ["Engineered-Cre", "Engineered-Flp"]
    group_donor = ["Donor-NRXN1del", "Donor-Ctrl"]
    group_types = [group_engineered, group_donor]

    df_info = pd.DataFrame(columns = ["time", "cell type", "condition target", "n_cells (target)", \
                                        "condition reference", "n_cells (reference)", "run_DE", "n_pos_DEGs", "n_neg_DEGs"])
    idx = 0
    # iteration for each kind of sample group
    for group in group_types:
        meta_sub_g = meta.loc[meta["sample_group"].isin(group), : ].copy()
        # iteration over time
        for time in np.unique(meta_sub_g[L_time]):
            meta_sub_t = meta_sub_g.loc[meta_sub_g[L_time] == time, : ]
            # iteration over cell type
            for cell_type in np.unique(meta_sub_t[L_cellT]):
                meta_sub_c = meta_sub_t.loc[meta_sub_t[L_cellT] == cell_type, : ]
                condition_dict = dict(meta_sub_c["sample_group"].value_counts())

                # number of cells in either group must be equal or greater than 5.
                if min(condition_dict[group[0]], condition_dict[group[1]]) > 4:
                    adata_sub = adata[list(meta_sub_c.index.values), : ].copy()              
                    sc.tl.rank_genes_groups(adata_sub, groupby = "sample_group", groups = [group[0]], reference = group[1], 
                                            method = method, n_genes = n_genes, pts = True, use_raw = False)
                    df_DEG = format_DEGs(adata_sub)
                    DEG_file_name = time + "_" + cell_type + "_" + group[0] + "_vs_" + group[1] + ".txt"
                    df_DEG.to_csv(DEG_file_name, sep="\t")

                    n_pos = df_DEG.loc[(df_DEG["logfoldchanges"] > 0) & (df_DEG["pvals_adj"] < 0.05), :].shape[0]
                    n_neg = df_DEG.loc[(df_DEG["logfoldchanges"] < 0) & (df_DEG["pvals_adj"] < 0.05), :].shape[0]

                    df_info.loc[idx] = [time, cell_type, group[0], group[1], condition_dict[group[0]], condition_dict[group[1]],\
                                                "Y", n_pos, n_neg]
                else:
                    df_info.loc[idx] = [time, cell_type, group[0], group[1], condition_dict[group[0]], condition_dict[group[1]],\
                                                "N", 0, 0]
                idx += 1
df_info.to_csv("cell_count_DEG_matrix.txt", sep = "\t")                

'''
    # 1. engineered
    adata_engineered = adata[adata.obs["sample_group"].isin(group_engineered), :].copy()
    print("engineered")
    for time in np.unique(adata_engineered.obs["timepoint"]):
        print(time)
        for cell_type in np.unique(adata_engineered.obs["cellclass_draft_v2"]):
            print(cell_type)
            adata_engineered_time_celltype = adata_engineered[(adata_engineered.obs["timepoint"] == time) & (adata_engineered.obs["cellclass_draft_v2"] == cell_type), :].copy()
            adata_engineered_time_celltype.obs["sample_group"].value_counts()
            try:
                print(adata_engineered_time_celltype.obs["sample_group"].value_counts())
                sc.tl.rank_genes_groups(adata_engineered_time_celltype, groupby = "sample_group", groups = ["Engineered-Cre"], reference = "Engineered-Flp", 
                                        method = method, n_genes = n_genes, pts = True, use_raw = False)
                df_DEG = format_DEGs(adata_engineered_time_celltype)
                DEG_file_name = time + "_" + cell_type + "_engineered_Flp_vs_Cre.txt"
                df_DEG.to_csv(DEG_file_name, sep="\t")
            except:
                print(cell_type)


    # 2. patient-derived
    adata_patient_derived = adata[adata.obs["sample_group"].isin(group_donor), :].copy()
    print("donor")
    for time in np.unique(adata_patient_derived.obs["timepoint"]):
        print(time)        
        for cell_type in np.unique(adata_patient_derived.obs["cellclass_draft_v2"]):
            print(cell_type)
            adata_patient_derived_time_celltype = adata_patient_derived[(adata_patient_derived.obs["timepoint"] == time) & (adata_patient_derived.obs["cellclass_draft_v2"] == cell_type), :].copy()
            print(adata_patient_derived_time_celltype.obs["sample_group"].value_counts())
            sc.tl.rank_genes_groups(adata_patient_derived_time_celltype, groupby = "sample_group", groups = ["Donor-NRXN1del"], reference = "Donor-Ctrl", 
                                    method = method, n_genes = n_genes, pts = True, use_raw = False)
            df_DEG = format_DEGs(adata_patient_derived_time_celltype)
            DEG_file_name = time + "_" + cell_type + "_patient_derived_NRXN1del_vs_ctrl.txt"
            df_DEG.to_csv(DEG_file_name, sep="\t")
'''