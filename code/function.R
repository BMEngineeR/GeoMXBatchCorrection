#' Calculate Silhoutte and kBET Score
#'
#' @param spe SpatialExperiment object.
#' @param batch_name Wanted or unwanted factor used to calculate. 
#'
#' @return A list with Silhoutte and kBET Score. 
#' @export
#'
#' @examples
cal_sil_kbet <- function(spe, batch_name = "patient", seed = 123){
  data <- assay(spe, i = 2)
  # find HVG
  dec <- scran::modelGeneVar(data)
  top_genes <- scran::getTopHVGs(dec, n = 1000)
  data_hvg <- t(data[rownames(data) %in% top_genes, ])
  
  # Silhoutte
  pca.data.f <- gmodels::fast.prcomp(data_hvg, center = TRUE)
  dd <- as.matrix(dist(pca.data.f$x[, 1:10]))
  batch <- colData(spe)[, batch_name]
  batch <- factor(
    batch,
    levels = sort(unique(batch)),
    labels = seq_along(sort(unique(batch)))
  )
  batch.silhouette <- summary(cluster::silhouette(as.numeric(batch), dd))$avg.width
  
  # kBET
  set.seed(seed)
  batch <- colData(spe)[, batch_name]
  k0 <- floor(mean(table(batch))) #neighbourhood size: mean batch size 
  knn <- get.knn(data_hvg, k = k0, algorithm = "cover_tree")
  batch.estimate <- kBET(data_hvg, batch, k = k0, knn = knn,plot = F)
  batch.kbet <- mean(batch.estimate$stats$kBET.observed)
  
  df_core <- data.frame(batch.silhouette, batch.kbet)
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
    .spe <- geomxNorm(spe, method = normalization_method)
  }
  .spe <- findNCGs(.spe, batch_name = unwanted_variable, top_n = ncg_topn)
  # [to do] change factors to biological variation remain based on meta data
  # here to setup wanted biology variable
  .spe <- geomxBatchCorrection(.spe, factors = wanted_variable,
                               NCGs = metadata(.spe)$NCGs, k = ruv4_k)
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
    dplyr::mutate(rank = rank(value_rank_scale))%>% 
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
      ggplot(aes(x = ruv4_k, y = value_rank_scale, color = factor(ncg_topn))) +
      geom_point(aes(shape = normalization_methods)) +
      geom_line(aes(group = paste(normalization_methods, ncg_topn))) +
      facet_grid(score ~ variable) +
      labs(title = "Individual Scaled Score") +
      scale_x_continuous(breaks = scales::breaks_width(1))
  }
  if (category == 1 & value == "rank") {
    p <- list_rank$rank %>%
      ggplot(aes(x = ruv4_k, y = rank, color = factor(ncg_topn))) +
      geom_point(aes(shape = normalization_methods)) +
      geom_line(aes(group = paste(normalization_methods, ncg_topn))) +
      facet_grid(score ~ variable) +
      labs(title = "Individual Score Ranking") +
      scale_x_continuous(breaks = scales::breaks_width(1)) 
  }
  if (category == 2 & value == "scaled") {
    p <- list_rank$mean_rank_variable %>%
      ggplot(aes(x = ruv4_k, y = value_rank_scale, color = factor(ncg_topn))) +
      geom_point(aes(shape = normalization_methods)) +
      geom_line(aes(group = paste(normalization_methods, ncg_topn))) +
      facet_wrap(~ variable, scales = "free") +
      labs(title = "Mean Scaled Score (per Variable)") +
      scale_x_continuous(breaks = scales::breaks_width(1)) 
  }
  if (category == 2 & value == "rank") {
    p <- list_rank$mean_rank_variable %>%
      ggplot(aes(x = ruv4_k, y = rank, color = factor(ncg_topn))) +
      geom_point(aes(shape = normalization_methods)) +
      geom_line(aes(group = paste(normalization_methods, ncg_topn))) +
      facet_wrap(~ variable, scales = "free") +
      labs(title = "Mean Scaled Score Ranking (per Variable)") +
      scale_x_continuous(breaks = scales::breaks_width(1))   
  }
  if (category == 3 & value == "scaled") {
    p <- list_rank$mean_rank_score %>%
      ggplot(aes(x = ruv4_k, y = value_rank_scale, color = factor(ncg_topn))) +
      geom_point(aes(shape = normalization_methods)) +
      geom_line(aes(group = paste(normalization_methods, ncg_topn))) +
      facet_wrap(~ score, scales = "free") +
      labs(title = "Mean Scaled Score (per Score)") +
      scale_x_continuous(breaks = scales::breaks_width(1))
  }
  if (category == 3 & value == "rank") {
    p <- list_rank$mean_rank_score %>%
      ggplot(aes(x = ruv4_k, y = rank, color = factor(ncg_topn))) +
      geom_point(aes(shape = normalization_methods)) +
      geom_line(aes(group = paste(normalization_methods, ncg_topn))) +
      facet_wrap(~ score, scales = "free") +
      labs(title = "Mean Scaled Score Ranking (per Score)") +
      scale_x_continuous(breaks = scales::breaks_width(1))
  }
  if (category == 4 & value == "scaled") {
    p <- list_rank$overall_rank %>%
      ggplot(aes(x = ruv4_k, y = value_rank_scale, color = factor(ncg_topn))) +
      geom_point(aes(shape = normalization_methods)) +
      geom_line(aes(group = paste(normalization_methods, ncg_topn))) + 
      labs(title = "Overall Mean Scaled Score") +
      scale_x_continuous(breaks = scales::breaks_width(1))          
  }         
  if (category == 4 & value == "rank") {
    p <- list_rank$overall_rank %>%
      ggplot(aes(x = ruv4_k, y = rank, color = factor(ncg_topn))) +
      geom_point(aes(shape = normalization_methods)) +
      geom_line(aes(group = paste(normalization_methods, ncg_topn))) + 
      labs(title = "Overall Mean Scaled Score Ranking") +
      scale_x_continuous(breaks = scales::breaks_width(1))           
  }
  return(p)
}

# Demo #####
if (FALSE) {
  raw_exp_mat <- read.csv("data/Masked_subsample_rawexp.csv", row.names = 1)
  meta_df <- read.csv("data/Masked_subsample_metadata.csv", row.names = 1)
  spe <- SpatialExperiment(assay = list(counts = raw_exp_mat), colData = meta_df)
  
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
}