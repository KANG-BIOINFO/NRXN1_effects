import scanpy as sc
import numpy as np
import pandas as pd
import pickle

adata = sc.read("combined_all_organoids_v2_doubletDetection.h5ad")
adata.obs.rename({"doublet":"doublet_doubletdetection","doublet_score":"doublet_score_doubletdetection"},axis = 1, inplace = True)
with open("doublet_scores_scrublet.pl", "rb") as f:
    prediction_scrublet = pickle.load(f)

adata.obs["doublet_scrublet"] = prediction_scrublet["predicted_doublets"]
adata.obs["doublet_score_scrublet"] = prediction_scrublet["doublet_scores"]
adata.write("combined_all_organoids_v2_doubletDetection_scrublet.h5ad")
adata.obs = adata.obs.astype({"doublet_scrublet": "category"})
sc.pl.umap(adata, color = ["doublet_doubletdetection", "doublet_score_doubletdetection", "doublet_scrublet", "doublet_score_scrublet"], save = "doublet.pdf")