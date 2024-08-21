#' Calculate Silhoutte and kBET Score
#'
#' @param spe SpatialExperiment object.
#' @param batch_name Wanted or unwanted factor used to calculate. 
#'
#' @return A list with Silhoutte and kBET Score. 
#' @export
#'
#' @examples
cal_sil_kbet <- function(spe, batch_name = "patient", seed = 123) {
  data <- SummarizedExperiment::assay(spe, i = 2)
  # find HVG
  dec <- scran::modelGeneVar(data)
  top_genes <- scran::getTopHVGs(dec, n = 1000)
  data_hvg <- t(data[rownames(data) %in% top_genes, ])
  
  # Silhoutte
  pca_data_f <- gmodels::fast.prcomp(data_hvg, center = TRUE)
  dd <- as.matrix(dist(pca_data_f$x[, 1:10]))
  batch <- SummarizedExperiment::colData(spe)[, batch_name]
  batch <- factor(
    batch,
    levels = sort(unique(batch)),
    labels = seq_along(sort(unique(batch)))
  )
  batch_silhouette <- summary(cluster::silhouette(as.numeric(batch), dd))$avg.width
  
  # kBET
  set.seed(seed)
  batch <- SummarizedExperiment::colData(spe)[, batch_name]
  k0 <- floor(mean(table(batch))) #neighbourhood size: mean batch size 
  knn <- FNN::get.knn(data_hvg, k = k0, algorithm = "cover_tree")
  batch_estimate <- kBET::kBET(data_hvg, batch, k = k0, knn = knn, plot = FALSE)
  batch_kbet <- mean(batch_estimate$stats$kBET.observed)
  
  df_core <- data.frame(batch_silhouette, batch_kbet)
  names(df_core) <- stringr::str_c(batch_name, c("silhouette", "kbet"), sep = "_")
  return(df_core)
}


#' Run GeoMx Batch Correction
#'
#' @param spe SpatialExperiment object. 
#' @param wanted_variable Wanted factor want to retain. 
#' @param unwanted_variable Unwanted factor want to remove. 
#' @param normalization_method Method used for normalization. 
#' @param ncg_topn Number of negative control gene used for batch correction. 
#' @param ruv4_k The k of RUV4. Based on RUV’s documentation, it is suggest to use the smallest k possible where the observed technical variation is no longer observed.
#'
#' @return A batch-corrected SpatialExperiment object. 
#' @export
#'
#' @examples
geomx_batch_correction <- function(spe, wanted_variable, unwanted_variable, 
                                   normalization_method, ncg_topn, ruv4_k) {
  if (normalization_method == "CPM") {
    .spe <- scater::logNormCounts(spe)
  } else {
    .spe <- standR::geomxNorm(spe, method = normalization_method)
  }
  .spe <- standR::findNCGs(.spe, batch_name = unwanted_variable, top_n = ncg_topn)
  # [to do] change factors to biological variation remain based on meta data
  # here to setup wanted biology variable
  .spe <- standR::geomxBatchCorrection(
    .spe, 
    factors = wanted_variable,
    NCGs = S4Vectors::metadata(.spe)$NCGs, 
    k = ruv4_k
  )
  return(.spe)
}


#' Score Batch Correction
#'
#' @param spe SpatialExperiment object. 
#' @param wanted_variable Wanted factor want to retain. 
#' @param unwanted_variable Unwanted factor want to remove. 
#' @param normalization_method Method used for normalization. 
#' @param ncg_topn Number of negative control gene used for batch correction. 
#' @param ruv4_k The k of RUV4. Based on RUV’s documentation, it is suggest to use the smallest k possible where the observed technical variation is no longer observed.
#'
#' @return A data.frame with Silhoutte and kBET scores. 
#' @export
#'
#' @examples
score_batch_correction <- function(spe, wanted_variable, unwanted_variable, 
                                   normalization_methods, ncg_topn, ruv4_k) {
  .spe <- geomx_batch_correction(spe, wanted_variable, unwanted_variable, 
                                 normalization_methods, ncg_topn, ruv4_k)
  .df_score <- c(unwanted_variable, wanted_variable) %>%
    purrr::map(~ cal_sil_kbet(.spe, batch_name = .x)) %>%
    dplyr::bind_cols()
  .df_parameter <- data.frame(
    normalization_methods = normalization_methods, 
    ncg_topn = ncg_topn, 
    ruv4_k = ruv4_k, 
    tag = stringr::str_c(normalization_methods, ncg_topn, ruv4_k, sep = "-")
  )
  df_score <- dplyr::bind_cols(.df_parameter, .df_score)
  return(df_score)
}
score_batch_correction_safely <- purrr::safely(score_batch_correction)


#' Rank Batch Correction
#'
#' @param df_score A data.frame with Silhoutte and kBET scores (generated from `score_batch_correction()` fucntion).
#' @param wanted_variable Wanted factor want to retain.
#'
#' @return A list with 4 data.frames (1) raw scaled score values and ranks, (2) mean scaled values of different score categories (per batch) and rank, (3) mean scaled values of different batched (per score) and rank, and (4) mean scaled values of different batches and score categories. 
#' @export
#'
#' @examples
rank_batch_correction <- function(df_score, wanted_variable) {
  # rank of scores, lower is better
  df_rank <- df_score %>% 
    dplyr::select(tag, dplyr::contains("_silhouette"), dplyr::contains("_kbet")) %>%
    tidyr::pivot_longer(-tag, names_to = "name", values_to = "value") %>%
    tidyr::separate(name, sep = "_(?=silhouette|kbet)", into = c("variable", "score")) %>% 
    # for wanted variables, larger is better; for unwanted variables, smaller is better
    dplyr::mutate(value_rank = ifelse(variable %in% wanted_variable, -value, value)) %>%
    dplyr::group_by(variable, score) %>% 
    dplyr::mutate(
      value_rank_scale = (value_rank - min(value_rank)) / (max(value_rank) - min(value_rank)), 
      # if min == max, the scaled value is set as 0.5
      value_rank_scale = ifelse(is.nan(value_rank_scale), 0.5, value_rank_scale), 
      rank = rank(value_rank_scale)
    ) %>%
    dplyr::ungroup()
  
  # mean rank of Silhouette and kBET score for each batch
  df_mean_rank_variable <- df_rank %>% 
    dplyr::group_by(tag, variable) %>% 
    dplyr::summarise(value_rank_scale = mean(value_rank_scale), .groups = "drop") %>% 
    dplyr::group_by(variable) %>%
    dplyr::mutate(rank = rank(value_rank_scale)) %>% 
    dplyr::ungroup()
  df_mean_rank_score <- df_rank %>% 
    dplyr::group_by(tag, score) %>% 
    dplyr::summarise(value_rank_scale = mean(value_rank_scale), .groups = "drop") %>% 
    dplyr::group_by(score) %>%
    dplyr::mutate(rank = rank(value_rank_scale)) %>% 
    dplyr::ungroup()
  
  # mean rank of each batch
  df_overall_rank <- df_rank %>%
    dplyr::group_by(tag) %>% 
    dplyr::summarise(value_rank_scale = mean(value_rank_scale), .groups = "drop") %>% 
    dplyr::mutate(rank = rank(value_rank_scale))
  
  list_rank <- list(
    rank = df_rank, 
    mean_rank_variable = df_mean_rank_variable, 
    mean_rank_score = df_mean_rank_score, 
    overall_rank = df_overall_rank
  ) 
  list_rank <- list_rank %>%
    purrr::map(
      ~ .x %>%
        tidyr::separate(tag, sep = "-", into = c("normalization_methods", "ncg_topn", "ruv4_k")) %>% 
        dplyr::mutate(dplyr::across(c(ncg_topn, ruv4_k), as.numeric))
    )
  return(list_rank)
}


#' Plot Performance of Batch Correction
#'
#' @param list_rank A list with 4 data.frames generated by `rank_batch_correction()` function 
#' @param category Categories to plot: 1 = raw scaled score values or ranks, 2 = mean scaled values of different score categories (per batch) or rank, 3 = mean scaled values of different batched (per score) or rank, and 4 = mean scaled values of different batches and score categories or rank.
#' @param value Value of y-axis: "scaled" = scaled score values, and "rank" = rank. 
#'
#' @return ggplot object. 
#' @export
#'
#' @examples
plot_batch_correction <- function(list_rank, category = 1, value = "scaled") {
  if (category == 1 & value == "scaled") {
    p <- list_rank$rank %>%
      ggplot2::ggplot(ggplot2::aes(x = ruv4_k, y = value_rank_scale, color = factor(ncg_topn))) +
      ggplot2::geom_point(ggplot2::aes(shape = normalization_methods)) +
      ggplot2::geom_line(ggplot2::aes(group = paste(normalization_methods, ncg_topn))) +
      ggplot2::facet_grid(score ~ variable) +
      ggplot2::labs(title = "Individual Scaled Score") +
      ggplot2::scale_x_continuous(breaks = scales::breaks_width(1))
  }
  if (category == 1 & value == "rank") {
    p <- list_rank$rank %>%
      ggplot2::ggplot(ggplot2::aes(x = ruv4_k, y = rank, color = factor(ncg_topn))) +
      ggplot2::geom_point(ggplot2::aes(shape = normalization_methods)) +
      ggplot2::geom_line(ggplot2::aes(group = paste(normalization_methods, ncg_topn))) +
      ggplot2::facet_grid(score ~ variable) +
      ggplot2::labs(title = "Individual Score Ranking") +
      ggplot2::scale_x_continuous(breaks = scales::breaks_width(1)) 
  }
  if (category == 2 & value == "scaled") {
    p <- list_rank$mean_rank_variable %>%
      ggplot2::ggplot(ggplot2::aes(x = ruv4_k, y = value_rank_scale, color = factor(ncg_topn))) +
      ggplot2::geom_point(ggplot2::aes(shape = normalization_methods)) +
      ggplot2::geom_line(ggplot2::aes(group = paste(normalization_methods, ncg_topn))) +
      ggplot2::facet_wrap(~ variable, scales = "free") +
      ggplot2::labs(title = "Mean Scaled Score (per Variable)") +
      ggplot2::scale_x_continuous(breaks = scales::breaks_width(1)) 
  }
  if (category == 2 & value == "rank") {
    p <- list_rank$mean_rank_variable %>%
      ggplot2::ggplot(ggplot2::aes(x = ruv4_k, y = rank, color = factor(ncg_topn))) +
      ggplot2::geom_point(ggplot2::aes(shape = normalization_methods)) +
      ggplot2::geom_line(ggplot2::aes(group = paste(normalization_methods, ncg_topn))) +
      ggplot2::facet_wrap(~ variable, scales = "free") +
      ggplot2::labs(title = "Mean Scaled Score Ranking (per Variable)") +
      ggplot2::scale_x_continuous(breaks = scales::breaks_width(1))   
  }
  if (category == 3 & value == "scaled") {
    p <- list_rank$mean_rank_score %>%
      ggplot2::ggplot(ggplot2::aes(x = ruv4_k, y = value_rank_scale, color = factor(ncg_topn))) +
      ggplot2::geom_point(ggplot2::aes(shape = normalization_methods)) +
      ggplot2::geom_line(ggplot2::aes(group = paste(normalization_methods, ncg_topn))) +
      ggplot2::facet_wrap(~ score, scales = "free") +
      ggplot2::labs(title = "Mean Scaled Score (per Score)") +
      ggplot2::scale_x_continuous(breaks = scales::breaks_width(1))
  }
  if (category == 3 & value == "rank") {
    p <- list_rank$mean_rank_score %>%
      ggplot2::ggplot(ggplot2::aes(x = ruv4_k, y = rank, color = factor(ncg_topn))) +
      ggplot2::geom_point(ggplot2::aes(shape = normalization_methods)) +
      ggplot2::geom_line(ggplot2::aes(group = paste(normalization_methods, ncg_topn))) +
      ggplot2::facet_wrap(~ score, scales = "free") +
      ggplot2::labs(title = "Mean Scaled Score Ranking (per Score)") +
      ggplot2::scale_x_continuous(breaks = scales::breaks_width(1))
  }
  if (category == 4 & value == "scaled") {
    p <- list_rank$overall_rank %>%
      ggplot2::ggplot(ggplot2::aes(x = ruv4_k, y = value_rank_scale, color = factor(ncg_topn))) +
      ggplot2::geom_point(ggplot2::aes(shape = normalization_methods)) +
      ggplot2::geom_line(ggplot2::aes(group = paste(normalization_methods, ncg_topn))) + 
      ggplot2::labs(title = "Overall Mean Scaled Score") +
      ggplot2::scale_x_continuous(breaks = scales::breaks_width(1))          
  }         
  if (category == 4 & value == "rank") {
    p <- list_rank$overall_rank %>%
      ggplot2::ggplot(ggplot2::aes(x = ruv4_k, y = rank, color = factor(ncg_topn))) +
      ggplot2::geom_point(ggplot2::aes(shape = normalization_methods)) +
      ggplot2::geom_line(ggplot2::aes(group = paste(normalization_methods, ncg_topn))) + 
      ggplot2::labs(title = "Overall Mean Scaled Score Ranking") +
      ggplot2::scale_x_continuous(breaks = scales::breaks_width(1))           
  }
  return(p)
}


#' Run PCA and UMAP for Before/After Batch Correction
#'
#' @param spe_b SpatialExperiment object before batch correction.
#' @param spe_a SpatialExperiment object after batch correction.
#' @param seed A number for random seed.
#'
#' @return A list of SpatialExperiment object of before/after batch correction with PCA and UMAP.
#' @export
#'
#' @examples
ba_pca_umap <- function(spe_b, spe_a, seed = 123) {
  set.seed(seed)
  if (!"logcounts" %in% SummarizedExperiment::assayNames(spe_b)) {
    spe_b <- scater::logNormCounts(spe_b)
  }
  spe_b <- scater::runPCA(spe_b)
  spe_a <- scater::runPCA(spe_a)
  spe_b <- scater::runUMAP(spe_b, dimred = "PCA")
  spe_a <- scater::runUMAP(spe_a, dimred = "PCA")
  ba_spe <- list(spe_b = spe_b, spe_a = spe_a)
  return(ba_spe)
}


#' Plot PCA and UMAP for Before/After Batch Correction
#'
#' @param ba_spe A list of SpatialExperiment object of before/after batch correction with PCA and UMAP (generated from `ba_pca_umap()`).
#' @param color_by Specification of a column metadata field or a feature to colour by.
#' @param shape_by Specification of a column metadata field or a feature to shape by.
#' @param point_size Size of points.
#'
#' @return A patchwork object of 4 figures.
#' @export
#'
#' @examples
ba_plot_pca_umap <- function(ba_spe, color_by, shape_by, point_size = 5) {
  spe_b <- ba_spe[["spe_b"]]
  spe_a <- ba_spe[["spe_a"]]
  p_1 <- scater::plotPCA(spe_b, colour_by = color_by, shape_by = shape_by, point_size = point_size) +
    ggplot2::labs(title = "Before Batch Correction") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  p_2 <- scater::plotUMAP(spe_b, colour_by = color_by, shape_by = shape_by, point_size = point_size)
  p_3 <- scater::plotPCA(spe_a, colour_by = color_by, shape_by = shape_by, point_size = point_size) +
    ggplot2::labs(title = "After Batch Correction") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  p_4 <- scater::plotUMAP(spe_a, colour_by = color_by, shape_by = shape_by, point_size = point_size)
  p <- patchwork::wrap_plots(p_1, p_2, p_3, p_4, byrow = FALSE, guides = "collect")
  return(p)
}


# Demo #####
if (FALSE) {
  library(magrittr)
  
  raw_exp_mat <- read.csv("data/Masked_subsample_rawexp.csv", row.names = 1)
  meta_df <- read.csv("data/Masked_subsample_metadata.csv", row.names = 1)
  spe <- SpatialExperiment::SpatialExperiment(assay = list(counts = raw_exp_mat), colData = meta_df)
  
  # define unwanted and wanted factors
  wanted_variable <- c("Disease", "celltype")
  unwanted_variable <- c("patientID")
  
  # method of normalization
  normalization_methods <- c("TMM", "CPM", "upperquartile", "sizefactor")
  # number of negative control gene
  ncg_topn <- seq(from = 100, to = 500, by = 200)
  # k of RUV4 (The number of unwanted factors to use. Can be 0. Required for the RUV4 method.)
  ruv4_k <- seq(from = 3, to = 10, by = 2)
  
  # Score Batch Correction
  df_parameter <- tidyr::expand_grid(normalization_methods, ncg_topn, ruv4_k) 
  list_score <- df_parameter %>% 
    purrr::pmap(
      ~ score_batch_correction_safely(spe, wanted_variable, unwanted_variable, ..1, ..2, ..3)$result, 
      .progress = TRUE
    )
  df_score <- dplyr::bind_rows(list_score)
  list_rank <- rank_batch_correction(df_score, wanted_variable)
  plot_batch_correction(list_rank, 1, "scaled")
  plot_batch_correction(list_rank, 2, "scaled")
  plot_batch_correction(list_rank, 3, "scaled")
  plot_batch_correction(list_rank, 4, "scaled")
  plot_batch_correction(list_rank, 1, "rank")
  plot_batch_correction(list_rank, 2, "rank")
  plot_batch_correction(list_rank, 3, "rank")
  plot_batch_correction(list_rank, 4, "rank")

  # dimensional reduction before and after batch correction
  wanted_variable <- c("Disease", "celltype")
  unwanted_variable <- c("patientID")
  normalization <- "upperquartile"
  ncg_topn <- 500
  ruv4_k <- 5
  spe_bc <- geomx_batch_correction(spe, wanted_variable, unwanted_variable, normalization, ncg_topn, ruv4_k)
  ba_spe <- ba_pca_umap(spe, spe_bc)
  ba_plot_pca_umap(ba_spe, color_by = "patientID", shape_by = "celltype", point_size = 3) +
    patchwork::plot_annotation(title = stringr::str_c(normalization, ncg_topn, ruv4_k, sep = ", "))
}