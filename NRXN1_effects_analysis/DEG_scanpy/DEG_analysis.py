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

    # 0. initialization
    adata_all = sc.read("/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/combined_all_organoids_v4.h5ad")

    adata_donor = adata_all[adata_all.obs["sample_group"].isin(["Donor-Ctrl","Donor-NRXN1del"])]
    adata_engineered = adata_all[adata_all.obs["sample_group"].isin(["Engineered-Cre","Engineered-Flp"])]


    method = "wilcoxon" # [‘logreg’, ‘t-test’, ‘wilcoxon’, ‘t-test_overestim_var’]
    groupby = "condition"

    # donor
    for cell_type in np.unique(adata_donor.obs["cellclass_draft_v2"]):
        adata_sub = adata_donor[adata_donor.obs["cellclass_draft_v2"] == cell_type, : ].copy()
        if adata_sub.shape[0] > 10:
            DEG_file_name = "./samplegroup_celltype_KO/DEG_donor_" + cell_type + "_NRXNdel_vs_ctrl.txt"
            sc.tl.rank_genes_groups(adata_sub, groupby = groupby, groups = ["NRXN_del"], reference = "Control",
                                    method = method, n_genes = adata_sub.shape[1], pts = True, use_raw = False)
            df_DEG = format_DEGs(adata_sub)
            df_DEG.to_csv(DEG_file_name, sep="\t")

    # engineered
    for cell_type in np.unique(adata_engineered.obs["cellclass_draft_v2"]):
        adata_sub = adata_engineered[adata_engineered.obs["cellclass_draft_v2"] == cell_type, : ].copy()
        if adata_sub.shape[0] > 10:
            DEG_file_name = "./samplegroup_celltype_KO/DEG_engineered_" + cell_type + "_NRXNdel_vs_ctrl.txt"
            sc.tl.rank_genes_groups(adata_sub, groupby = groupby, groups = ["NRXN_del"], reference = "Control",
                                    method = method, n_genes = adata_sub.shape[1], pts = True, use_raw = False)
            df_DEG = format_DEGs(adata_sub)
            df_DEG.to_csv(DEG_file_name, sep="\t")