Raw_exp_mat <- read.csv("data/Masked_subsample_rawexp.csv", row.names = 1)
meta_df <- read.csv("data/Masked_subsample_metadata.csv", row.names = 1)


#' Calculate Silhoutte and kBET Score
#'
#' @param spe SpatialExperiment object.
#' @param batch.name Wanted or unwanted factor used to calculate. 
#'
#' @return A list with Silhoutte and kBET Score. 
#' @export
#'
#' @examples
cal_sil_kbet <- function(spe, batch.name = "patient"){
  data <- assay(spe, i = 2)
  # find HVG
  dec <- scran::modelGeneVar(data)
  top_genes <- scran::getTopHVGs(dec, n = 1000)
  data_hvg <- t(data[rownames(data) %in% top_genes, ])
  
  # Silhoutte
  pca.data.f <- gmodels::fast.prcomp(data_hvg, center = TRUE)
  dd <- as.matrix(dist(pca.data.f$x[, 1:10]))
  batch <- colData(spe)[, batch.name]
  batch <- factor(
    batch,
    levels = sort(unique(batch)),
    labels = seq_along(sort(unique(batch)))
  )
  batch.silhouette <- summary(cluster::silhouette(as.numeric(batch), dd))$avg.width
  
  # kBET
  batch <- colData(spe)[, batch.name]
  k0 <- floor(mean(table(batch))) #neighbourhood size: mean batch size 
  knn <- get.knn(data_hvg, k = k0, algorithm = "cover_tree")
  batch.estimate <- kBET(data_hvg, batch, k = k0, knn = knn,plot = F)
  batch.kbet <- mean(batch.estimate$stats$kBET.observed)
  
  df_core <- data.frame(batch.silhouette, batch.kbet)
  names(df_core) <- stringr::str_c(batch.name, c("silhouette", "kbet"), sep = "_")
  return(df_core)
}



#' Run GeoMx Batch Correction
#'
#' @param spe SpatialExperiment object. 
#' @param WantedVariable Wanted factor want to retain. 
#' @param unwantedVariable Unwanted factor want to remove. 
#' @param normalization_method Method used for normalization. 
#' @param findNCG_topn Number of negative control gene used for batch correction. 
#' @param BatchCorrection_k The k of RUV4. Based on RUV’s documentation, it is suggest to use the smallest k possible where the observed technical variation is no longer observed.
#'
#' @return A batch-corrected SpatialExperiment object. 
#' @export
#'
#' @examples
run_geomxBatchCorrection <- function(spe, WantedVariable, unwantedVariable, 
                                     normalization_method, findNCG_topn, BatchCorrection_k) {
  if (normalization_method == "CPM") {
    .spe <- scater::logNormCounts(spe)
  } else {
    .spe <- geomxNorm(spe, method = normalization_method)
  }
  .spe <- findNCGs(.spe, batch_name = unwantedVariable, top_n = findNCG_topn)
  # [to do] change factors to biological variation remain based on meta data
  # here to setup wanted biology variable
  .spe <- geomxBatchCorrection(.spe, factors = WantedVariable,
                               NCGs = metadata(.spe)$NCGs, k = BatchCorrection_k)
  return(.spe)
}

evaluate_BatchCorrection <- function(spe, WantedVariable, unwantedVariable, 
                                     normalization_methods, findNCG_topn, BatchCorrection_k) {
  .spe <- run_geomxBatchCorrection(spe, WantedVariable, unwantedVariable, 
                                   normalization_methods, findNCG_topn, BatchCorrection_k)
  .df_score <- c(unwantedVariable, WantedVariable) %>%
    purrr::map(~ cal_sil_kbet(.spe, batch.name = .x)) %>%
    dplyr::bind_cols()
  .df_parameter <- data.frame(
    normalization_methods = normalization_methods, 
    findNCG_topn = findNCG_topn, 
    BatchCorrection_k = BatchCorrection_k, 
    tag = stringr::str_c(normalization_methods, findNCG_topn, BatchCorrection_k, sep = "-")
  )
  df_score <- dplyr::bind_cols(.df_parameter, .df_score)
  return(df_score)
}
evaluate_BatchCorrection_safely <- purrr::safely(evaluate_BatchCorrection)



# create spatial Experiment obj
spe <- SpatialExperiment(assay = list(counts = Raw_exp_mat), colData = meta_df)

# define unwanted factor
unwantedVariable <- c("patientID")
# define wanted factor
WantedVariable <- c("Disease","celltype")
normalization_methods <- c("TMM", "CPM", "upperquartile", "sizefactor")
# number of negative control gene
findNCG_topn <- seq(from = 100, to = 2000, by = 200)
# k of RUV4 (The number of unwanted factors to use. Can be 0. This is required for the RUV4 method.)
BatchCorrection_k <- seq(from = 3, to = 20, by = 2)

df_parameter <- expand.grid(
  normalization_methods = normalization_methods, 
  findNCG_topn = findNCG_topn, 
  BatchCorrection_k = BatchCorrection_k, 
  # score_batch = c(unwantedVariable, WantedVariable), 
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(tag = stringr::str_c(normalization_methods, findNCG_topn, BatchCorrection_k, sep = "-"))
set.seed(199681)
score_list <- vector("list", nrow(df_parameter))
progress <- txtProgressBar(min = 0, max = nrow(df_parameter), style = 3)
for (i in seq_along(score_list)) {
  .df_parameter <- df_parameter[i, ]
  normalization_methods <- .df_parameter$normalization_methods
  findNCG_topn <- .df_parameter$findNCG_topn
  BatchCorrection_k <- .df_parameter$BatchCorrection_k
  score_list[[i]] <- evaluate_BatchCorrection_safely(spe, WantedVariable, unwantedVariable, 
                                  normalization_methods, findNCG_topn, BatchCorrection_k)$result
  setTxtProgressBar(progress, i)
}
close(progress)
df_score <- dplyr::bind_rows(score_list)

df_score %>% readr::write_rds("20240818_BatchCorrection.rds")

df_score %>% 
  dplyr::select(tag, dplyr::contains("_silhouette"), dplyr::contains("_kbet")) %>%
  tidyr::pivot_longer(-tag, names_to = "name", values_to = "value") %>%
  tidyr::separate(name, sep = "_", into = c("variable", "score")) %>% 
  dplyr::mutate(value = ifelse(variable %in% WantedVariable, -value, value)) %>%
  dplyr::group_by(variable, score) %>% 
  dplyr::mutate(rank = rank(value)) %>% 
  dplyr::group_by(tag, variable) %>% 
  dplyr::summarise(meanrank = mean(rank)) %>%
  dplyr::group_by(tag) %>% 
  dplyr::summarise(rank = mean(meanrank)) %>% 
  dplyr::arrange(rank)



dplyr::bind_cols(df_parameter, df_score) %>% 
  tidyr::pivot_longer(c(dplyr::contains("_silhouette"), dplyr::contains("_kbet")), names_to = "name", values_to = "value") %>%
  tidyr::separate(name, sep = "_", into = c("variable", "score")) %>% 
  dplyr::mutate(value = ifelse(variable %in% WantedVariable, -value, value)) %>% 
  ggplot(aes(x = BatchCorrection_k, y = value, color = paste0(normalization_methods, findNCG_topn))) +
  geom_point() +
  geom_line() +
  facet_wrap(score ~ variable, scales = "free")
