import scanpy as sc
import numpy as np 
import pandas as pd 
import scrublet as scr
import pickle

adata = sc.read("../combined_all_organoids_v2.h5ad")
scrub = scr.Scrublet(adata.raw.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
combined = {"doublet_scores": doublet_scores,
            "predicted_doublets": predicted_doublets}

with open("doublet_scores.pl", "wb") as f:
    pickle.dump(combined, f)
