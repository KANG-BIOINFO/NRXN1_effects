import pandas as pd 
import numpy as np
import os 
import matplotlib.pyplot as plt 

outFolder_time = "time_samplegroup_KO"
outFolder_time_celltype = "time_samplegroup_celltype_KO"

df_demo = pd.read_csv("time_samplegroup_KO/d50_engineered_Cre_vs_Flp.txt", sep = "\t", header = 0, index_col = 0)

df_deg_values = pd.DataFrame(index = list(df_demo.index.values))
df_deg_logits = pd.DataFrame(index = list(df_demo.index.values))
df_colname_metadata = pd.DataFrame(columns = ["sample_group_type", "timepoint", "cell_type"])

for outFolder in [outFolder_time, outFolder_time_celltype]:
    files_de = [(outFolder + "/" + i) for i in os.listdir(outFolder) if (i.startswith("d") and (i.endswith(".txt")))]
    for file_de in files_de:
        # load data
        df_deg = pd.read_csv(file_de, sep = "\t", header = 0, index_col = 0).loc[list(df_demo.index.values), : ]
        colname = file_de.split("/")[1].rstrip(".txt")
        print(colname)

        # fill in values of adjusted p values
        df_deg_values[colname] = df_deg["pvals_adj"]

        # fill in values of up/down regulated genes
        df_deg_logits[colname] = 0
        up_genes = list(df_deg.loc[(df_deg["logfoldchanges"] > 0) & (df_deg["pvals_adj"] < 0.05), :].index.values)
        down_genes = list(df_deg.loc[(df_deg["logfoldchanges"] < 0) & (df_deg["pvals_adj"] < 0.05), :].index.values)
        df_deg_logits.loc[up_genes, colname] = 1
        df_deg_logits.loc[down_genes, colname] = -1

        # fill in column metadata
        sample_group = "Engineered" if "Cre" in colname else "Patient-derived"
        time = colname.split("_")[0]

        if outFolder == outFolder_time:
            df_colname_metadata.loc[colname] = [sample_group, time, ""]
        else:
            df_colname_metadata.loc[colname] = [sample_group, time, colname.split("_")[1]]

df_deg_values.to_csv("survey/DEG_values.txt", sep = "\t")
df_deg_logits.to_csv("survey/DEG_directions.txt", sep = "\t")
df_colname_metadata.to_csv("survey/DEG_colname_metadata.txt", sep = "\t")



