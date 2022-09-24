library(Seurat)
library(glmnet)
library(ggplot2)
library(cowplot)
library(foreach)
library(doMC)

#
getPopulationOffset = function(y){
  if(!is.factor(y))
    y=factor(y)
  if(length(levels(y))!=2)
    stop("y must be a two-level factor")
  off = sum(y==levels(y)[2])/length(y)
  off = log(off/(1-off))
  return(rep(off,length(y)))
}


#' Do the OvR fit for every variable.  This just does a simple CV selection of regularisation amount.  Far from ideal, but should be good enough for the main conclusions.
multinomialFitCV = function(x,y,nParallel=1,...){
  fits = list()
  if(nParallel>1)
    registerDoMC(cores=nParallel)
  #Do them in order of size
  marks = names(sort(table(as.character(y))))
  for(mark in marks){
    message(sprintf("Fitting model for variable %s",mark))
    fac = factor(y==mark)
    #The two main modes of failure are too few positives and errors constructing lambda.  These should be handled semi-gracefully
    fits[[mark]] = tryCatch(
      cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',intercept=FALSE,alpha=0.99,nfolds=10,type.measure='class',parallel=nParallel>1,...),
      error = function(e) {
        tryCatch(
          cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',intercept=FALSE,alpha=0.99,nfolds=10,type.measure='class',parallel=nParallel>1,lambda=exp(seq(-10,-3,length.out=100)),...),
          error = function(e) {
            warning(sprintf("Could not fit model for variable %s",mark))
            return(NULL)
          })
      })
  }
  return(fits)
}

# 2. train the model
train_predict_fromdata = function(query, ref_data_file, name) # query_data and data_file are all seurat data
{
  # load data
  query_dat = t(query@assays$RNA@data)
  
  ref = readRDS(ref_data_file)
  ref = NormalizeData(ref)
  ref_dat = t(ref@assays$RNA@data)
  
  # overlap gene
  overlap_genes = intersect(colnames(query_dat), colnames(ref_dat))
  query_dat = query_dat[, overlap_genes]
  ref_dat = ref_dat[, overlap_genes]
  
  # train the model
  classes = ref@meta.data[["Cell.Type"]]
  if (name == "Zhong"){classes = ref@meta.data[["cell_types"]]}
  fits = multinomialFitCV(ref_dat, classes, nParallel = 10) # the model
  
  # fit the model
  preds = list()
  for (mark in names(fits)){
    message(sprintf("Predicting probabilities for cluster %s",mark))
    preds[[mark]] = predict(fits[[mark]], 
                            newx = query_dat[, overlap_genes],
                            s = 'lambda.1se', newoffset = rep(0,nrow(query_dat)))
  }
  df_pred = do.call(cbind, preds)
  colnames(df_pred) = names(preds)
  df_pred = as.data.frame(df_pred)

  write.table(df_pred, paste0("prediction_output/output_cell_", name, ".txt"), sep = "\t")
  
  # produce a final output
  clusters = query$cellclass_draft_v2
  output = list()
  for (cluster in unique(clusters))
  {
    cells = names(clusters)[clusters == cluster]
    df_pred_sub = df_pred[cells,]
    output[[cluster]] = colMeans(df_pred_sub)
  }
  df_output = as.data.frame(do.call(cbind, output))
  write.table(df_output, paste0("prediction_output/output_cluster_", name, ".txt"), sep = "\t")
  
  # return the model
  return(fits)
}

# run the function
input_files = list.files("/data/aronow/Kang/single_cell_projects/10X-Pak/reference/raw_data/extract_subset/",
                         pattern = ".rds")

organoid = readRDS("../combined_all_organoids_filtered.rds")
df_cells = read.table("../ctype_v4.txt", sep = "\t", header = 1, row.names = 1)
df_cells = df_cells[colnames(organoid),]
organoid$cellclass_draft_v2 = df_cells$cellclass_draft_v2

# four organoid datasets
names = c("Kanton", "Paulsen", "Tanaka","Velasco")
models = list()
for (idx in 1:4)
{
  path = paste0("/data/aronow/Kang/single_cell_projects/10X-Pak/reference/raw_data/extract_subset/", input_files[idx])
  models[[idx]] = train_predict_fromdata(organoid, path, names[idx])
}

# one fetal brain data
name = "Zhong"
path = "/data/aronow/Kang/single_cell_projects/10X-Pak/reference/fetalBrain_Zhong_nature/FetalBrain_Zhong_raw_annotated.rds"
models[[5]] = train_predict_fromdata(organoid, path, name)





