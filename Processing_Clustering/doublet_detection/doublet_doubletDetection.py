import scanpy as sc
import numpy as np 
import pandas as pd 
import doubletdetection
import pickle

sc.settings.n_jobs=8
adata = sc.read("../combined_all_organoids_v2.h5ad")

'''
scrub = scr.Scrublet(adata.raw.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
combined = {"doublet_scores": doublet_scores,
            "predicted_doublets": predicted_doublets}
'''

clf = doubletdetection.BoostClassifier(n_iters = 25, use_phenograph=False, standard_scaling=True)
# raw_counts is a cells by genes count matrix
doublets = clf.fit(adata.raw.X).predict(p_thresh=1e-16, voter_thresh=0.5)
# higher means more likely to be doublet
doublet_score = clf.doublet_score()

adata.obs["doublet"] = doublets
adata.obs["doublet_score"] = doublet_score
adata.write("combined_all_organoids_v2_doubletDetection.h5ad")

with open("doublet_scores_doubletdetection.pl", "wb") as f:
    pickle.dump(scores, f)
