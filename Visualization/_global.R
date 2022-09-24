# genes
genes_mature = c('STMN2', 'GAP43', 'DCX')
genes_progenitor = c('VIM', 'HES1', 'SOX2')
genes_cycle = c('CCNB2', 'CDK1', 'MKI67')

genes_npc = c('PTTG1', 'PCNA', 'HMGB2', 'MKI67')
genes_choroid = c('TPBG')
genes_mesenchymal = c('DCN')
genes_microglia = c('CD74')
genes_glia = c('HOPX', 'FAM107A', 'TNC', 'AQP4', 'SLC1A3', 'PMP2')
genes_ipc = c('EOMES', 'NHLH1')
genes_opc = c('HES6', 'OLIG1', 'OLIG2')
genes_in = c('DLX1', 'DLX2', 'DLX5', 'SLC32A1', 'GAD1', 'GAD2')
genes_cn = c('TBR1', 'NEUROD2', 'NEUROD6', 'IRX3', 'FOXP2')

genes_atlas = c('VIM', 'HES1', 'PTTG1','MKI67', 'HOPX', 'AQP4', 'EOMES', 'HES6',
                'DLX1', 'DLX2', 'GAD1', 'NEUROD2', 'NEUROD6', 'FOXP2', 
                'TPBG', 'DCN', 'CD74')

# Cell metadata
cellclass_ordered = c("Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                      "Astroglia", "Glia cell", "oRG", "Intermediate cells", "Choroid Plexus",
                      "Mesenchymal cell", "Microglia",
                      "OPC", "IPC1","IPC2","IPC3",
                      "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                      "CN1","CN2","CN3","CN4","CN5") # no 'low sequencing depth' and 'Unknown'

cellclass_ordered_v2 = c("Cycling Dorsal NEC", "Cycling Ventral NEC", "Cycling NEC","Non-cycling NEC",
                         "Astroglia", "Glia cell", "oRG", "Intermediate cells", "Choroid Plexus",
                         "OPC", "IPC1","IPC2","IPC3",
                         "IN1","IN2","IN3","IN4","IN5","IN6","IN7",
                         "CN1","CN2","CN3","CN4","CN5") # no 'Mesenchymal cell' and 'Microglia'

colour_cellclass = c("#00008B","#0000CD","#4169E1","#4682B4",
                     "#79BA69","#E5E350","#5C8F47","#81C6C5","#919145",
                     "#97B68A","#708747",
                     "#F6CB41", "#C8AF85","#EEAF58","#D8BB6C",
                     "#B26FAA","#DA478E","#A36CA9","#BA408D","#A8376F","#E3CCCB","#DE759B",
                     "#CD4D2C","#E58A41","#D93D2F","#81231B","#C34D46")

colour_cellclass_v2 = c("#00008B","#0000CD","#4169E1","#4682B4",
                         "#79BA69","#E5E350","#5C8F47","#81C6C5","#919145",
                         "#F6CB41", "#C8AF85","#EEAF58","#D8BB6C",
                         "#B26FAA","#DA478E","#A36CA9","#BA408D","#A8376F","#E3CCCB","#DE759B",
                         "#CD4D2C","#E58A41","#D93D2F","#81231B","#C34D46")

sample_group_ordered = c('Engineered-Flp-d023', 'Engineered-Cre-d023',
                         'Engineered-Flp-d050', 'Engineered-Cre-d050',
                         'Engineered-Flp-d112', 'Engineered-Cre-d112',
                         'Donor-Ctrl-d022', 'Donor-NRXN1del-d022',
                         'Donor-Ctrl-d050', 'Donor-NRXN1del-d050',
                         'Donor-Ctrl-d101', 'Donor-NRXN1del-d101')
sample_group_ordered_v2 = c('Donor-Ctrl-d022', 'Donor-NRXN1del-d022',
                             'Donor-Ctrl-d050', 'Donor-NRXN1del-d050',
                             'Donor-Ctrl-d101', 'Donor-NRXN1del-d101',
                            'Engineered-Flp-d112', 'Engineered-Cre-d112')
  
colour_timepoint = c('#E5FFCC', '#66B2FF', '#FF6666')
colour_timepoint_v2 = c('#006600', '#3399FF', '#CC0000')

colour_genotype = c('#99CCFF', '#FF6666')
