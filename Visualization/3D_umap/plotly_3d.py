import plotly.express as px
import pandas as pd
import numpy as np

########## donor #############
df = pd.read_csv('./umap_3d_metadata_donor.txt', sep = '\t', header = 0, index_col = 0)
df.loc[df['cellclass_draft_v2'] == 'OPC', 'cellclass_draft_v2'] = 'IPC4'

# order cell type and extract important cells
cellclass_ordered_v2 = ["Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                         "Astroglia", "Glia cell", "oRG", "Intermediate cells", "Choroid Plexus",
                         "IPC4", "IPC1","IPC2","IPC3",
                         "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                         "CN1","CN2","CN3","CN4","CN5"] # no 'Mesenchymal cell' and 'Microglia'
colour_cellclass_v2 = ["#00008B","#0000CD","#4169E1","#4682B4",
                         "#79BA69","#E5E350","#5C8F47","#81C6C5","#919145",
                         "#F6CB41", "#C8AF85","#EEAF58","#D8BB6C",
                         "#B26FAA","#DA478E","#A36CA9","#BA408D","#A8376F","#E3CCCB","#DE759B",
                         "#CD4D2C","#E58A41","#D93D2F","#81231B","#C34D46"]
map_ = dict(zip(cellclass_ordered_v2, colour_cellclass_v2))
df = df.loc[df['cellclass_draft_v2'].isin(cellclass_ordered_v2), :].copy()

# remove other cells in 2d umap
df = df.loc[df['UMAP_1'] < 7.5, :].copy()
df['color'] = [map_[i] for i in df['cellclass_draft_v2']]
df['size'] = 0.4

fig = px.scatter_3d(df, x='umap_3d_1', y='umap_3d_2', z='umap_3d_3',
                    size = 'size',
                    size_max = 7.5,
                    color='cellclass_draft_v2',
                    color_discrete_map = map_
                   )
fig.write_html('v2_donor/donor_3d.html',include_plotlyjs='cdn')


########## all ###########
df = pd.read_csv('./umap_3d_metadata_all.txt', sep = '\t', header = 0, index_col = 0)
df.loc[df['cellclass_draft_v2'] == 'OPC', 'cellclass_draft_v2'] = 'IPC4'

# order cell type and extract important cells
cellclass_ordered_v2 = ["Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                         "Astroglia", "Glia cell", "oRG", "Intermediate cells", "Choroid Plexus",
                         "IPC4", "IPC1","IPC2","IPC3",
                         "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                         "CN1","CN2","CN3","CN4","CN5"] # no 'Mesenchymal cell' and 'Microglia'
colour_cellclass_v2 = ["#00008B","#0000CD","#4169E1","#4682B4",
                         "#79BA69","#E5E350","#5C8F47","#81C6C5","#919145",
                         "#F6CB41", "#C8AF85","#EEAF58","#D8BB6C",
                         "#B26FAA","#DA478E","#A36CA9","#BA408D","#A8376F","#E3CCCB","#DE759B",
                         "#CD4D2C","#E58A41","#D93D2F","#81231B","#C34D46"]
map_ = dict(zip(cellclass_ordered_v2, colour_cellclass_v2))
df = df.loc[df['cellclass_draft_v2'].isin(cellclass_ordered_v2), :].copy()

df['color'] = [map_[i] for i in df['cellclass_draft_v2']]
df['size'] = 0.4

fig = px.scatter_3d(df, x='umap_3d_1', y='umap_3d_2', z='umap_3d_3',
                    size = 'size',
                    size_max = 7.5,
                    color='cellclass_draft_v2',
                    color_discrete_map = map_
                   )
fig.write_html('v2_donor/all_3d.html',include_plotlyjs='cdn')


########## engineered ###########
df = pd.read_csv('./umap_3d_metadata_engineered.txt', sep = '\t', header = 0, index_col = 0)
df.loc[df['cellclass_draft_v2'] == 'OPC', 'cellclass_draft_v2'] = 'IPC4'

# order cell type and extract important cells
cellclass_ordered_v2 = ["Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                         "Astroglia", "Glia cell", "oRG", "Intermediate cells", "Choroid Plexus",
                         "IPC4", "IPC1","IPC2","IPC3",
                         "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                         "CN1","CN2","CN3","CN4","CN5"] # no 'Mesenchymal cell' and 'Microglia'
colour_cellclass_v2 = ["#00008B","#0000CD","#4169E1","#4682B4",
                         "#79BA69","#E5E350","#5C8F47","#81C6C5","#919145",
                         "#F6CB41", "#C8AF85","#EEAF58","#D8BB6C",
                         "#B26FAA","#DA478E","#A36CA9","#BA408D","#A8376F","#E3CCCB","#DE759B",
                         "#CD4D2C","#E58A41","#D93D2F","#81231B","#C34D46"]
map_ = dict(zip(cellclass_ordered_v2, colour_cellclass_v2))
df = df.loc[df['cellclass_draft_v2'].isin(cellclass_ordered_v2), :].copy()

df['color'] = [map_[i] for i in df['cellclass_draft_v2']]
df['size'] = 0.4

fig = px.scatter_3d(df, x='umap_3d_1', y='umap_3d_2', z='umap_3d_3',
                    size = 'size',
                    size_max = 7.5,
                    color='cellclass_draft_v2',
                    color_discrete_map = map_
                   )
fig.write_html('v2_donor/engineered_3d.html',include_plotlyjs='cdn')