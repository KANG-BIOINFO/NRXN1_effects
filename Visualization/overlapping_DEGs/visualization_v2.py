import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

df_num = pd.read_csv('num_pairwise_intersection_genes_donor.txt',sep='\t',header=0,index_col=0)
df_score = pd.read_csv('hypergeometric_logfdr_pairwise_intersection_genes_donor_v2.txt',sep='\t',header=0,index_col=0)
df_meta_all = pd.read_csv('num_pairwise_intersection_genes_donor_metadata_v2.txt',sep='\t',header=0,index_col=0)

cell_classes = ["Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                "Astroglia", "Glia cell", "oRG", "Intermediate cells",
                "OPC", "IPC1","IPC2","IPC3",
                "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                "CN1","CN2","CN3","CN4","CN5"]

colour_cellclass_v2 = ["#00008B","#0000CD","#4169E1","#4682B4",
                         "#79BA69","#E5E350","#5C8F47","#81C6C5",
                         "#F6CB41", "#C8AF85","#EEAF58","#D8BB6C",
                         "#B26FAA","#DA478E","#A36CA9","#BA408D","#A8376F","#E3CCCB","#DE759B",
                         "#CD4D2C","#E58A41","#D93D2F","#81231B","#C34D46"]
                         
dot_cmap = dict(zip(cell_classes, colour_cellclass_v2))
Norm = mpl.colors.Normalize(vmin = 0, vmax = 100, clip = True)

timepoints = ['d22','d50','d101']
timepoints_v2 = ['d022','d050','d101']

for timepoint, timepoint_v2 in zip(timepoints, timepoints_v2):
    
    df_meta = df_meta_all.loc[df_meta_all['timepoint'] == timepoint_v2, : ].copy()
    print(df_meta.shape)

    cell_classes_avail = [i for i in cell_classes if i in df_meta['cell_class'].unique()]
    items_avail = ['Donor-NRXN1del_' + timepoint + '_' + i for i in cell_classes_avail]

    # df_score = df_score.loc[items_avail, : ].copy() # reorder rows and columns
    # df_score = df_score[items_avail].copy()

    n_classes = len(cell_classes_avail)

    print(cell_classes_avail)
    fig, ax = plt.subplots()

    # for i in range(n_classes):
    #     for j in range(0, n_classes - i):
    #         item_name_i = 'Donor-NRXN1del_' + timepoint + '_' + cell_classes_avail[i]
    #         item_name_j = 'Donor-NRXN1del_' + timepoint + '_' + cell_classes_avail[i+j]
    #         num = df_num.loc[item_name_i, item_name_j]
    #         score = df_score.loc[item_name_i, item_name_j]
            
    #         ax.scatter(np.array(i), np.array(n_classes - j - i - 1), s = num / 5, c = score, norm = Norm, cmap = 'Blues')

    xs = []; ys = []; nums = []; scores = []
    for i in range(n_classes):
        for j in range(i, n_classes):
            # print('i: ',i, ', j: ', j)
            item_name_i = 'Donor-NRXN1del_' + timepoint + '_' + cell_classes_avail[i]
            item_name_j = 'Donor-NRXN1del_' + timepoint + '_' + cell_classes_avail[j]

            # print(item_name_i, item_name_j)
            num = df_num.loc[item_name_i, item_name_j]
            score = df_score.loc[item_name_i, item_name_j]
            
            xs.append(i)
            ys.append(n_classes - j - 1)
            nums.append(num)
            scores.append(score)

    dots = ax.scatter(xs, ys, s = np.array(nums) / 5, c = scores, norm = Norm, cmap = 'Blues')
    plt.colorbar(dots, location = 'right', orientation = 'vertical', shrink = 0.3, label = '-log10(FDR)')

    # draw lines
    for i in range(n_classes):
        j = n_classes - i
        ax.plot([i-0.5, i-0.5], [-0.5, j-0.5], linewidth = 0.6, color = 'grey', alpha = 0.6)
        ax.plot([-0.5, j-0.5], [i-0.5, i-0.5], linewidth = 0.6, color = 'grey', alpha = 0.6)

    # add axis symbols
    for i in range(n_classes):
        ax.scatter(-1, n_classes -i - 1, s = 25, color = dot_cmap[cell_classes_avail[i]], label = cell_classes_avail[i]) # y axis symbols
        ax.scatter(i, n_classes -i, s = 25, color = dot_cmap[cell_classes_avail[i]]) # x axis symbols (top)

    legend1 = plt.legend(*dots.legend_elements('sizes', num = 5), bbox_to_anchor = (0.7, 1))
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, l, bbox_to_anchor = (1, 1), prop = {'size':5})
    
    plt.gca().add_artist(legend1)

    ax.axis('off')
    fig.savefig('donor_' + timepoint + '_overlapDEGs_v2.pdf', bbox_inches = 'tight')