library(standR)
library(SpatialExperiment)
library(limma)
library(scater)
library(ExperimentHub)
library(ggalluvial)
require(FNN)
library(kBET)
# setup file
setwd("/bmbl_data/yuzhou/collaborative/Sizun_lab/EBV_project/Raw Data/batch effect correction/")

# read exp
# this should be expression count file (gene by ROI), the count matrix should be after QC (followed by standR or GeomxTools).
# DSP-1001660020565-D-C06 DSP-1001660020565-D-C07 DSP-1001660020565-D-C08 DSP-1001660020565-D-C09 DSP-1001660020565-D-C10
# A2M                        38                      68                      41                       6                      21
# NAT2                        6                      38                      10                       5                       8
# ACADM                      15                      49                       9                       9                       7
# ACADS                      15                      66                      25                      12                      17
# ACAT1                       3                      29                      11                       7                       4
Raw_exp_mat <- read.csv("../ebv_chl_rawcounts_filtered.csv",row.names = 1)
colnames(Raw_exp_mat) <- gsub("[.]","-",colnames(Raw_exp_mat))
# read meta (example show as below)
# TMA EBV Subtype Sex   Age Source   Source2               Sample_ID  segment   patient
# DSP-1001660020565-D-C06 DFCI neg      NS   M Young   Mass left neck DSP-1001660020565-D-C06    Tumor H528-1neg
# DSP-1001660020565-D-C07 DFCI neg      NS   M Young   Mass left neck DSP-1001660020565-D-C07   CD4mem H528-1neg
# DSP-1001660020565-D-C08 DFCI neg      NS   M Young   Mass left neck DSP-1001660020565-D-C08 CD4naive H528-1neg
# DSP-1001660020565-D-C09 DFCI neg      NS   M Young   Mass left neck DSP-1001660020565-D-C09   CD8mem H528-1neg
# DSP-1001660020565-D-C10 DFCI neg      NS   M Young   Mass left neck DSP-1001660020565-D-C10 CD8naive H528-1neg
meta_df <- read.csv("MetaData_CT.csv",row.names = 1)
# check meta data ROI name is as same as raw count ROI name
identical(colnames(Raw_exp_mat),rownames(meta_df))
# create spatial Experiment obj
spe <- SpatialExperiment(assay = list(counts = Raw_exp_mat), colData = meta_df)
# generate diff normalization + batch combination
# here you can specify the parameter for grid search
normalization_methods <- c("TMM", "CPM", "upperquartile", "sizefactor")
findNCG_topn <- c(100,300,500,600,700,900,1100,1300,1500)
BatchCorrection_k <- c(3,5,7,9,11)
length(BatchCorrection_k)*length(findNCG_topn)*length(normalization_methods)
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
    tmp_spe <- findNCGs(tmp_spe, batch_name = "patient", top_n = tmp_topn)
    
    for (k in 1:length(BatchCorrection_k)) {
      tmp_k <- BatchCorrection_k[k]
      # [to do] change factors to biological variation remain 
      tmp_spe <- geomxBatchCorrection(tmp_spe, factors = c("segment", "EBV"),
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

# save object list as checkpoint
library(qs)
# qsave(spe.list, file = "spe.list.qs")
spe.list <- qread("spe.list.qs")
# Silhoutte
cal_sil <- function(se.object = spe.list[[names(spe.list)[1]]],
                    obj.name = names(spe.list)[1],
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
                     obj.name = names(spe.list)[1],
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

# run for each of factor, including patient (batch confounder), segment (biological variation), and EBV (biological variation)
batch.name = "patient"
sil_score_list <- list()
k_bet_list <- list()
# progress bar
total_iterations <- length(spe.list)

# Initialize progress bar
progress <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration_count <- 0
for (i in 1:length(spe.list)){
  se.object = spe.list[[names(spe.list)[i]]]
  obj.name = names(spe.list)[i]
  sil_score <- cal_sil(se.object = se.object, obj.name = obj.name ,batch.name = batch.name )
  k_bet.obj <- cal_kbet(se.object = se.object, obj.name = obj.name ,batch.name = batch.name )
  sil_score_list <- c(sil_score_list, sil_score)
  k_bet_list <- c(k_bet_list, mean(k_bet.obj$stats$kBET.observed))
  names(sil_score_list)[length(sil_score_list)] <- names(k_bet_list)[length(k_bet_list)] <-  paste0(batch.name
                                                                                                    ,"_",obj.name)
  # Update progress bar
  iteration_count <- iteration_count + 1
  setTxtProgressBar(progress, iteration_count)
}
# save resulst and name it
qsave( k_bet_list ,file = paste0(batch.name,"_k_bet_list.qs"))
qsave( sil_score_list ,file = paste0(batch.name,"_sil_score_list.qs"))

#################################################################################################
# evaluation of all factors
# minimize
patient_kbet <- unlist(qread("patient_k_bet_list.qs"))
patient_sil <- unlist(qread("patient_sil_score_list.qs"))
identical(names(patient_sil),names(patient_kbet))
names(patient_sil) <- names(patient_kbet) <- gsub("^patient_","",names(patient_kbet))

# max 
segment_kbet <- unlist(qread("segment_k_bet_list.qs"))
segment_sil <- unlist(qread("segment_sil_score_list.qs"))
identical(names(segment_kbet),names(segment_sil))
names(segment_kbet) <- names(segment_sil) <- gsub("^segment_","",names(segment_sil))
segment_kbet <- segment_kbet[names(patient_sil)]
segment_sil <- segment_sil[names(patient_sil)]
# max 
EBV_kbet <- unlist(qread("EBV_k_bet_list.qs"))
EBV_sil <- unlist(qread("EBV_sil_score_list.qs"))
identical(names(EBV_kbet),names(EBV_sil))
names(EBV_kbet) <- names(EBV_sil) <- gsub("^patient_","",names(EBV_kbet))
identical(names(patient_kbet),names(segment_kbet))
# increasing order of the confounder (the lower the better)
patient_kbet_rank <- rank(patient_kbet)
patient_sil_rank <- rank(patient_sil)
# mean rank of batch
mean_rank_batch <- (patient_kbet_rank + patient_sil_rank)/2
sort(mean_rank_batch)
# decreasing order of the biological variance (the higher the better)
segment_kbet_rank <- rank(-segment_kbet)
segment_sil_rank <- rank(-segment_sil)
#
mean_rank_bio_CT <- (segment_kbet_rank + segment_sil_rank)/2
# decreasing order of biological variance (the higher the better)
EBV_kbet_rank <- rank(-EBV_kbet)
EBV_sil_rank <- rank(-EBV_sil)
#
mean_rank_bio_EBV <- (EBV_kbet_rank + EBV_sil_rank)/2

# consider bio and batch
score_table <- data.frame(row.names = gsub("^patient_","",names(mean_rank_batch)),
                          patient_kbet = patient_kbet,
                          patient_kbet_rank = patient_kbet_rank,
                          patient_sil = patient_sil,
                          patient_sil_rank = patient_sil_rank,
                          CT_kbet = segment_kbet,
                          CT_kbet_rank = segment_kbet_rank,
                          CT_sil = segment_sil,
                          CT_sil_rank = segment_sil_rank,
                          EBV_sil = EBV_sil,
                          EBV_sil_rank = EBV_sil_rank,
                          EBV_kbet = EBV_kbet,
                          EBV_kbet_rank = EBV_kbet_rank,
                          mean_rank_patient = mean_rank_batch,
                          mean_rank_CT = mean_rank_bio_CT,
                          mean_rank_EBV  = mean_rank_bio_EBV,
                          overall = (mean_rank_batch+ mean_rank_bio_CT+mean_rank_bio_EBV)/3)
score_table_sub <- score_table[,c("mean_rank_batch",
                                  "mean_rank_bio_CT",
                                  "mean_rank_bio_EBV",
                                  "mean_rank_overall")]
write.csv(score_table, file = "score_table_grid_search.csv",quote = F)
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



