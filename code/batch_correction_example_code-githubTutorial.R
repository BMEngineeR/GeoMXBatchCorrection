library(standR)
library(SpatialExperiment)
library(limma)
library(scater)
library(ExperimentHub)
library(ggalluvial)
require(FNN)
library(kBET)
# read exp
# this should be expression count file (gene by ROI), the count matrix should be after QC (followed by standR or GeomxTools).
# DSP-1001660020565-D-C06 DSP-1001660020565-D-C07 DSP-1001660020565-D-C08 DSP-1001660020565-D-C09 DSP-1001660020565-D-C10
# A2M                        38                      68                      41                       6                      21
# NAT2                        6                      38                      10                       5                       8
# ACADM                      15                      49                       9                       9                       7
# ACADS                      15                      66                      25                      12                      17
# ACAT1                       3                      29                      11                       7                       4
Raw_exp_mat <- read.csv("/bmbl_data/yuzhou/collaborative/Sizun_lab/EBV_project/Code/GitHubBatchEffectRemoval/GeoMXBatchCorrection/data/Masked_subsample_rawexp.csv",row.names = 1)
# read meta (example show as below)
# TMA EBV Subtype Sex   Age Source   Source2               Sample_ID  segment   patient
# DSP-1001660020565-D-C06 DFCI neg      NS   M Young   Mass left neck DSP-1001660020565-D-C06    Tumor H528-1neg
# DSP-1001660020565-D-C07 DFCI neg      NS   M Young   Mass left neck DSP-1001660020565-D-C07   CD4mem H528-1neg
# DSP-1001660020565-D-C08 DFCI neg      NS   M Young   Mass left neck DSP-1001660020565-D-C08 CD4naive H528-1neg
# DSP-1001660020565-D-C09 DFCI neg      NS   M Young   Mass left neck DSP-1001660020565-D-C09   CD8mem H528-1neg
# DSP-1001660020565-D-C10 DFCI neg      NS   M Young   Mass left neck DSP-1001660020565-D-C10 CD8naive H528-1neg
meta_df <- read.csv("/bmbl_data/yuzhou/collaborative/Sizun_lab/EBV_project/Code/GitHubBatchEffectRemoval/GeoMXBatchCorrection/data/Masked_subsample_metadata.csv",row.names = 1)

# create spatial Experiment obj
spe <- SpatialExperiment(assay = list(counts = Raw_exp_mat), colData = meta_df)
# generate diff normalization + batch combination
# here you can specify the parameter for grid search
unwantedVariable <- c("patientID")
WantedVariable <- c("Disease","celltype")
normalization_methods <- c("TMM", "CPM", "upperquartile", "sizefactor")
findNCG_topn <- seq(from = 100, to = 100, by = 200)
BatchCorrection_k <- seq(from = 3, to = 3, by = 2)
print(paste0("Combination Parameter Number: ",length(BatchCorrection_k)*length(findNCG_topn)*length(normalization_methods)))
# RUV4 parameters
# factors: the factor of interest, i.e. the biological variation to keep;
# NCGs: the list of negative control genes detected using the function findNCGs;
# k: the number of unwanted factors to use. Based on RUV’s documentation, it is suggest to use the smallest k possible where the observed technical variation is no longer observed.
spe.list <- list()
# Total number of iterations
total_iterations <- length(normalization_methods) * length(findNCG_topn) * length(BatchCorrection_k)

# Initialize progress bar
progress <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration_count <- 0

for (n_m in 1:length(normalization_methods)) {
  tmp_n_m <- normalization_methods[n_m]
  tmp_spe <- geomxNorm(spe, method = tmp_n_m)
  if(tmp_n_m == "CPM"){tmp_spe <- scater::logNormCounts(spe)}
  
  for (topn in 1:length(findNCG_topn)) {
    tmp_topn <- findNCG_topn[topn]
    #### here need to set up the batch factor based on meta data
    tmp_spe <- findNCGs(tmp_spe, batch_name = unwantedVariable, top_n = tmp_topn)
    
    for (k in 1:length(BatchCorrection_k)) {
      tmp_k <- BatchCorrection_k[k]
      # [to do] change factors to biological variation remain based on meta data
      # here to setup wanted biology variable
      tmp_spe <- geomxBatchCorrection(tmp_spe, factors = WantedVariable,
                                      NCGs = metadata(tmp_spe)$NCGs, k = k)
      tmp_name <- paste0(tmp_n_m, "-nNCG_", tmp_topn, "-k_", tmp_k)
      spe.list <- c(spe.list, tmp_spe)
      names(spe.list)[length(spe.list)] <- tmp_name
      
      # Update progress bar
      iteration_count <- iteration_count + 1
      setTxtProgressBar(progress, iteration_count)
    }
  }
}

# Close progress bar
close(progress)

# Silhoutte
cal_sil <- function(se.object = spe.list[[names(spe.list)[1]]],
                    #obj.name = names(spe.list)[1],
                    batch.name = "patient"){
  data <- assay(se.object,i = 2)
  # find HVG
  dec <-scran:: modelGeneVar(data)
  top_genes <- scran::getTopHVGs(dec, n = 1000)
  data_hvg <- t(data[rownames(data) %in% top_genes,])
  #
  batch <- colData(se.object)[,batch.name]
  pca.data.f <- gmodels::fast.prcomp(data_hvg, center=TRUE)
  dd <- as.matrix(dist(pca.data.f$x[, 1:10]))
  batch.silhouette <- summary(cluster::silhouette(as.numeric(factor(batch,
                                                                    levels = sort(unique(batch)),
                                                                    labels = 1:length(sort(unique(batch))))), dd))$avg.width
  return(batch.silhouette)
}

# kBET
cal_kbet <- function(se.object = spe.list[[names(spe.list)[1]]],
                    # obj.name = names(spe.list)[1],
                     batch.name = "patient"){
  data <- assay(se.object,i = 2)
  # find HVG
  dec <-scran:: modelGeneVar(data)
  top_genes <- scran::getTopHVGs(dec, n = 1000)
  data_hvg <- t(data[rownames(data) %in% top_genes,])
  #
  batch <- colData(se.object)[,batch.name]
  k0=floor(mean(table(batch))) #neighbourhood size: mean batch size 
  knn <- get.knn(data_hvg, k=k0, algorithm = 'cover_tree')
  batch.estimate <- kBET(data_hvg, batch, k = k0, knn = knn,plot = F)
  return(batch.estimate)
}

# run for each of factor, including patientID (batch confounder), celltype (biological variation), and Disease (biological variation)

batch.name = c(unwantedVariable, WantedVariable)
sil_score_list <- list()
k_bet_list <- list()
# progress bar
total_iterations <- length(spe.list)*length(batch.name)

# Initialize progress bar
progress <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration_count <- 0
for (i in 1:length(spe.list)){
  se.object = spe.list[[names(spe.list)[i]]]
  obj.name = names(spe.list)[i]
  for (j in 1:length(batch.name)){
    sil_score <- cal_sil(se.object = se.object, batch.name = batch.name[j] )
    k_bet.obj <- cal_kbet(se.object = se.object, batch.name = batch.name[j] )
    sil_score_list <- c(sil_score_list, sil_score)
    k_bet_list <- c(k_bet_list, mean(k_bet.obj$stats$kBET.observed))
    names(sil_score_list)[length(sil_score_list)] <-  names(k_bet_list)[length(k_bet_list)] <- paste0(batch.name[j],"_",obj.name)
    # Update progress bar
    iteration_count <- iteration_count + 1
    setTxtProgressBar(progress, iteration_count)
  }
}


#################################################################################################
# Xiaojie: you need to wrap up this function a bit (e.g., ranking all factor including unwanted variable and wanted variable)

# evaluation of all factors
# minimize unwanted variance
unwantedVariable_df <- matrix(rep(NA,length(unwantedVariable)*length(spe.list)*5),ncol =length(unwantedVariable)*5 )

unwantedVariable_combinations <- expand.grid(unwantedVariable,  c("_sil","_kbet","_Silrank","_kbetrank","_MeanRank"))

colnames(unwantedVariable_df) <- paste0(unwantedVariable_combinations$Var1, unwantedVariable_combinations$Var2)

rownames(unwantedVariable_df) <-  names(spe.list)
for (i in 1:length(unwantedVariable)){
  tmp_unwantedVariable <- unwantedVariable[i]
  tmp_kbet_name_vec <- grep(tmp_unwantedVariable,(names(k_bet_list)),value = F)
  tmp_kbet <-  unlist(k_bet_list[tmp_kbet_name_vec])
  tmp_sil_name_vec <- grep(tmp_unwantedVariable,(names(sil_score_list)),value = F)
  tmp_sil <-  unlist(sil_score_list[tmp_sil_name_vec])
  
  #
  names(tmp_sil) <- names(tmp_kbet) <- gsub(paste0("^",tmp_unwantedVariable,"_"),"",names(tmp_sil))
  tmp_unwantedVariable_kbet_rank <- rank(tmp_kbet)
  tmp_unwantedVariable_sil_rank <- rank(tmp_sil)
  # mean rank of batch
  tmp_unwantedVariable_kbet_rank <- tmp_unwantedVariable_kbet_rank[rownames(unwantedVariable_df)]
  tmp_unwantedVariable_sil_rank <- tmp_unwantedVariable_sil_rank[rownames(unwantedVariable_df)]
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_sil"))] <- as.numeric(tmp_sil)
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_kbet"))] <- as.numeric(tmp_kbet)
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_Silrank"))] <- as.numeric(tmp_unwantedVariable_sil_rank)
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_kbetrank"))] <- as.numeric(tmp_unwantedVariable_kbet_rank)
  # cal mean rank of sil and kbet
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_MeanRank"))] <- rowSums(unwantedVariable_df[,c(paste0(tmp_unwantedVariable,c("_Silrank")),
                                                                                                paste0(tmp_unwantedVariable,c("_kbetrank")))])/2
}

# Maximize wanted variance
WantedVariable_df <- matrix(rep(NA,length(WantedVariable)*length(spe.list)*5),ncol =length(WantedVariable)*5 )
WantedVariable_combinations <- expand.grid(WantedVariable, c("_sil","_kbet","_Silrank","_kbetrank","_MeanRank"))

colnames(WantedVariable_df) <- paste0(WantedVariable_combinations$Var1, WantedVariable_combinations$Var2)
rownames(WantedVariable_df) <-  names(spe.list)
for (i in 1:length(WantedVariable)){
  tmp_WantedVariable <- WantedVariable[i]
  tmp_kbet_name_vec <- grep(tmp_WantedVariable,(names(k_bet_list)),value = F)
  tmp_kbet <-  unlist(k_bet_list[tmp_kbet_name_vec])
  tmp_sil_name_vec <- grep(tmp_WantedVariable,(names(sil_score_list)),value = F)
  tmp_sil <-  unlist(sil_score_list[tmp_sil_name_vec])
  
  #
  names(tmp_sil) <- names(tmp_kbet) <- gsub(paste0("^",tmp_WantedVariable,"_"),"",names(tmp_sil))
  tmp_WantedVariable_kbet_rank <- rank(-tmp_kbet)
  tmp_WantedVariable_sil_rank <- rank(-tmp_sil)
  # mean rank of batch
  tmp_WantedVariable_kbet_rank <- tmp_WantedVariable_kbet_rank[rownames(WantedVariable_df)]
  tmp_WantedVariable_sil_rank <- tmp_WantedVariable_sil_rank[rownames(WantedVariable_df)]
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_sil"))] <- as.numeric(tmp_sil)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_kbet"))] <- as.numeric(tmp_kbet)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_Silrank"))] <- as.numeric(tmp_WantedVariable_sil_rank)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_kbetrank"))] <- as.numeric(tmp_WantedVariable_kbet_rank)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_MeanRank"))] <- rowSums(WantedVariable_df[,c(paste0(tmp_WantedVariable,c("_Silrank")),
                                                                                               paste0(tmp_WantedVariable,c("_kbetrank")))])/2
}

# merge results of wanted variable and unwanted variable

merge_res <- cbind.data.frame(unwantedVariable_df,WantedVariable_df)
merge_res <- merge_res[,grep("_MeanRank",colnames(merge_res))]
merge_res$overall_rank <- rowSums(merge_res[,grep("_MeanRank",colnames(merge_res))])/length(merge_res[,grep("_MeanRank",colnames(merge_res))]) 
merge_res <- merge_res[order(merge_res$overall_rank),]
##################################################################################################
# validate the best parameter
# c("TMM", "CPM", "upperquartile", "sizefactor")
# tmp_spe <- geomxNorm(spe, method = "upperquartile")
set.seed(123)
tmp_spe <- scater::logNormCounts(spe)
tmp_spe <- findNCGs(tmp_spe, batch_name = "patient", top_n = 1500)

tmp_spe <- geomxBatchCorrection(tmp_spe, factors = c("segment", "EBV"),
                                NCGs = metadata(tmp_spe)$NCGs, k = 11)
dec <-scran:: modelGeneVar(tmp_spe)
top_genes <- scran::getTopHVGs(dec, n = 1000)
exp_mat <- assay(tmp_spe,i = 2)

library(Seurat)
Seurat_obj <- CreateSeuratObject(exp_mat[top_genes,])
Seurat_obj <- AddMetaData(object = Seurat_obj,metadata = as.data.frame(colData(tmp_spe)))
Seurat_obj  <- FindVariableFeatures(Seurat_obj,nfeatures = 1000)
Seurat_obj  <- ScaleData(Seurat_obj)
Seurat_obj <- RunPCA(Seurat_obj)
Seurat_obj <- RunUMAP(Seurat_obj,dims = 1:20)
p1 <- DimPlot(Seurat_obj,group.by = "patient",pt.size = 3)+ theme(legend.position="none")
p2 <- DimPlot(Seurat_obj,group.by = "segment",pt.size = 3)
p3 <- DimPlot(Seurat_obj,group.by = "EBV",pt.size = 3)
p1+p2+p3
# write the post batch correction table
write.csv(exp_mat,file = "EBV_chl_batchcorrected-CPM-nNCGH_1500-K_11.csv",quote = F)



