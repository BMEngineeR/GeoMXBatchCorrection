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

# subsmaple data
table(meta_df$patient)
subsample_ID <- c("H528-10pos", "H528-11neg", "H528-12neg", "H528-12pos", "H528-13neg", "H528-14neg")
meta_df_sub <- meta_df[meta_df$patient %in% subsample_ID, c("Sample_ID","EBV","patient","segment")]
dim(meta_df_sub)
meta_df_sub[1:5,]
Raw_exp_mat_sub <- Raw_exp_mat[,meta_df_sub$Sample_ID]
Raw_exp_mat_sub[1:5,1:5]
colnames(Raw_exp_mat_sub) <- paste0("sample",1:ncol(Raw_exp_mat_sub))
rownames(meta_df_sub) <- meta_df_sub$Sample_ID <- colnames(Raw_exp_mat_sub)
meta_df_sub[1:5,]
meta_df_sub$Disease <- meta_df_sub$EBV
meta_df_sub$EBV <- NULL 
meta_df_sub[1:5,]
meta_df_sub$celltype <- meta_df_sub$segment
meta_df_sub$segment <- NULL
meta_df_sub[1:5,]

meta_df_sub$patientID <- factor(meta_df_sub$patient,levels = unique(meta_df_sub$patient), labels = c(paste0("Patient",1:6)))
table(meta_df_sub$patientID)
meta_df_sub$patient <- NULL
meta_df_sub
write.csv(meta_df_sub,"Masked_subsample_metadata.csv",sep = ",",quote = F)
write.csv(Raw_exp_mat_sub,"Masked_subsample_rawexp.csv",sep = ",",quote = F)



