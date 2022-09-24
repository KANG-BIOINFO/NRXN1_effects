import numpy as np 
import pandas as pd 
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt 
import seaborn as sns

# matplotlib.rcParams['text.usetex']=True
# matplotlib.rcParams['text.latex.unicode']=True

# load information
adata = sc.read('/data/aronow/Kang/single_cell_projects/10X-Pak/integration/patient_derived_engineered/v2_filtered/combined_all_organoids_v5.h5ad')
markers = ['VIM', 'HES1', 'PTTG1','MKI67', 'HOPX', 'AQP4', 'EOMES', 'HES6',
            'DLX1', 'DLX2', 'GAD1', 'NEUROD2', 'NEUROD6', 'FOXP2', 
            'TPBG', 'DCN', 'CD74']

markers_v2 = ['DHFR', 'TPX2', 'ASPM', 'CDC20','PCNA', 
             'NTRK2', 'EDNRB', 'WLS', 'TPBG','NHLH1', 
             'DCX', 'NCAM1', 'STMN4','SLC32A1', 'GAD1', 'GAD2',
             'DLX1', 'DLX2', 'DLX5', 'TBR1', 'NEUROD2', 'NEUROD6']

# proposed by ChangHui
markers_v3 = ['S100B', 'SLC1A3', 'AQP4', 'GFAP', 'SOX9', 'SPARCL1', # ASTROCYTES
            'ALDOC', 'HES5', 'CA2', 'AGT', 'AQP4', 'IL33', 'HEPN1', 'WNT7B', # GLIA
            'TNC', 'HOPX', 'MOXD1', 'FAM107A', 'PTPRZ1', 'LIFR', 'STAT3', 'PAX6', 'CLU', 'CARHSP1', 'SEMA5A', # OUTER RG
            'HES1', 'CRYAB', 'FBXO32', 'PALLD', #VENTRICULAR RG
            'PAX6','SOX2','PDGFD','GLI3','VIM','HES1',#RG
            'OLIG1','OLIG2','SOX10',#OPC
            'EOMES','NEUROG1','PPP1R17','BTG2','NEUROD1','PENK',#IPC
            'SLC17A7','SLC17A6','EMX1', #CN
            'CUX2','INHBA','BTG1','LPL','CITED2','PLXND1','SATB2','TBR1','BHLHE22','SYBU','LHX2','CUX1','POU3F3','LMO4','ZEB2',#CALLOSAL PROJECTION NEURONS
            'FEZF2','CNTN6','BCL11B','CRYM','LDB2','SOX5','TLE4','LMO3', # CORTICOSPINAL MOTOR NEURONS 
            'LHX6','NKX2-1','ARX','SOX6','SATB1','NPAS1','NR2F2','OTX2','ETV1','LHX6','LHX8','SOX6',#INTERNEURON
            'GAD1','GAD2','CALB1','SST','VIP','CALB2','NPY','RELN','CCK']#INTERNEURON

markers_v4 = ['DHFR', 'TPX2', 'ASPM', 'CDC20','PCNA', 'S100B', 'FAM107A', 'WNT7B', 'SOX9',
             'NTRK2', 'EDNRB', 'SLC1A3', 'WLS', 'HES1', 'ETV1', 'TPBG', 'OTX2', 'SOX6',  
             'NEUROG1', 'NHLH1', 'NEUROD1',
             'DCX', 'NCAM1', 'STMN4', 'SLC17A6','SLC32A1', 'GAD1', 'GAD2',
             'DLX1', 'DLX2', 'DLX5', 'ARX', 'TBR1', 'NEUROD2', 'NEUROD6',
             'LDB2','BCL11B','FEZF2','LHX2','SATB2','LPL','SLC17A7']

cellclass_ordered_v2 = ["Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                         "Astroglia", "Glia cell", "oRG", "Choroid Plexus", "Intermediate cells",
                         "OPC", "IPC1","IPC2","IPC3",
                         "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                         "CN1","CN2","CN3","CN4","CN5"]
# subsetting                        
'''adata = adata[adata.obs['cellclass_draft_v2'].isin(cellclass_ordered_v2), : ].copy()
adata.obs['cellclass_draft_v2'].cat.reorder_categories(cellclass_ordered_v2, inplace = True)
print(len(np.intersect1d(adata.var_names, markers)))
print(len(markers))
sc.pl.stacked_violin(adata, var_names = markers_v2, groupby = 'cellclass_draft_v2', 
                    use_raw = False, rotation = 90, row_palette = sns.color_palette("pastel"),
                    dendrogram = False, save = 'test_stackviolin_scanpy_2-5.pdf')'''

adata = adata[adata.obs['cellclass_draft_v2'].isin(cellclass_ordered_v2), : ].copy()
adata.obs['cellclass_draft_v2'].cat.reorder_categories(cellclass_ordered_v2, inplace = True)
print(len(np.intersect1d(adata.var_names, markers)))
print(len(markers))
sc.pl.stacked_violin(adata, var_names = markers_v3, groupby = 'cellclass_draft_v2', 
                    use_raw = False, rotation = 90, row_palette = sns.color_palette("pastel"),
                    dendrogram = False, save = 'test_stackviolin_scanpy_2-6.pdf')

sc.pl.stacked_violin(adata, var_names = markers_v4, groupby = 'cellclass_draft_v2', 
                    use_raw = False, rotation = 90, row_palette = sns.color_palette("pastel"),
                    dendrogram = False, save = 'test_stackviolin_scanpy_2-8.pdf')