# Utility functions for Saunders, Srivatsan, et al. Nature (2023)

# extract coordinates for a PAGA graph from a cell data set 
get_paga_graph <- function(cds, reduction_method = "UMAP") {
  
  cluster_result <- cds@clusters[[reduction_method]]$cluster_result
  
  cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g,
                                                     cluster_result$optim_res,
                                                     qval_thresh = 0.05, 
                                                     FALSE)
  
  cluster_g <- cluster_graph_res$cluster_g
  cluster_g <- igraph::set_vertex_attr(cluster_g, "name", value = stringr::str_replace(igraph::V(cluster_g)$name, "cell_membership", ""))
  cluster_g
}


# add 2d or 3d umap coordinates to your coldata
append_umap_coordinates = function(cds, umap_3D = F){
  
  if (!umap_3D){
    colData(cds)$umap1 = reducedDim(x = cds,
                                    type = "UMAP")[,1]
    colData(cds)$umap2 = reducedDim(x = cds,
                                    type = "UMAP")[,2]
  }
  if (umap_3D){
    colData(cds)$umap3d_1 = reducedDim(x = cds,
                                       type = "UMAP")[,1]
    colData(cds)$umap3d_2 = reducedDim(x = cds,
                                       type = "UMAP")[,2]
    colData(cds)$umap3d_3 = reducedDim(x = cds,
                                       type = "UMAP")[,3]
  }
  return(cds)
}

# calculate hotspot function
calc_hotspot <- function(cds, 
                         compare_col, compare,
                         subset_col = NULL, subset = NULL, 
                         reduction_method="PCA") {
  
  if (!(is.null(subset_col) & is.null(subset))) {
    cds = cds[, replace_na(colData(cds)[[subset_col]] == subset, F)]
  }
  pData(cds)$Cell <- pData(cds)$cell
  df = pData(cds) %>% as.data.frame()
  spatial_weights = calculateSpatialWeights(cds, reduction_method = reduction_method)
  lg <- calculateLocalG(cds = cds, lw = spatial_weights$lw, 
                        wc = spatial_weights$wc, column_name = compare_col, var = compare)
  lg.df <- mergeLocalG(lg, df, compare_col) %>% getPval("bonferroni", "twosided")
  lg.df
}

# input gene short names and get associated gene ids
get.gene.ids = function (cds, names) {
  return(as.character(rowData(cds)[rowData(cds)$gene_short_name %in% names, "id"]))
}

#' calculate spatial weights
#' functions take from graph_test
#' @k is how many neighbors in knn graph
#' @reduction_method default is UMAP
#'
calculateSpatialWeights <- function(cds, k=15, reduction_method = "UMAP"){
  
  lw <- monocle3:::calculateLW(cds = cds, k = k,
                               verbose = FALSE,
                               neighbor_graph = "knn",
                               reduction_method = reduction_method,
                               list(method = "nn2"))
  
  wc <- spdep::spweights.constants(lw, zero.policy = TRUE, adjust.n = TRUE)
  
  return(list(lw=lw, wc=wc))
}


#' @var is a specific label value
#' @column_name of labels of interest
#' @lw is output of calculateLW
calculateLocalG <- function(cds, lw, wc, column_name, var) {
  
  # absence or presence of specified variable
  # df$LG <- ifelse(df[,column_name] == var, 1,0)
  # z = df[,"LG"]
  # names(z) <- df$Cell
  
  cds$LG <- ifelse(pData(cds)[,column_name] == var, 1,0)
  z = pData(cds)[,"LG"]
  names(z) <- pData(cds)$Cell
  
  # calculate local g
  localG = spdep::localG(x=z, listw = lw, wc$n, wc$S0, alternative = NULL, zero.policy = TRUE, spChk = F)
  localG.df = localG[1:length(localG)] %>% as.data.frame()
  colnames(localG.df)[1] <- as.character(var)
  rownames(localG.df) = pData(cds)$Cell
  localG.df
}

#' select local g values of interest + merge with full df
#' only values 1 mean anything in the local G calculation
mergeLocalG <- function(localG_result, df, column) {
  
  localG_result$Cell = row.names(localG_result)
  wt = localG_result %>% merge(df, by=c("Cell"))
  wt = wt[wt[,column] == colnames(localG_result)[1],]
  colnames(wt)[2] <- "localG"
  wt
}

#' to turn the zscore into a pvalue
#' method : c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
getPval <- function(df, method, tail) {
  
  # I think this is right, but should check it...
  if (tail=="onesided"){
    df = df %>% mutate(pval = pnorm(-abs(localG)))
  } else if (tail=="twosided"){
    df = df %>% mutate(pval = 2*pnorm(-abs(localG)))
  }
  adjusted_pval = p.adjust( p = df$pval, method = method)
  df %>% mutate("adjusted_pval" = adjusted_pval)
}


#' compare col : column in cds that is your categorical label, ex: "gene_target"
#' compare: which category you want to determine hotspot, ex: "tbxta"
#' subset col: if you want to do hotspot by broad grouping, ex: "cell_type_broad"
#' subset: the specific subset, ex: "fast muscle"
calc_hotspot <- function(cds,
                         compare_col, compare,
                         subset_col = NULL, subset = NULL,
                         reduction_method="UMAP",
                         method = "bonferroni",
                         tail = "twosided") {
  
  if (!(is.null(subset_col) & is.null(subset))) {
    cds = cds[, tidyr::replace_na(colData(cds)[[subset_col]] == subset, F)]
  }
  pData(cds)$Cell <- pData(cds)$cell
  df = pData(cds) %>% as.data.frame()
  spatial_weights = calculateSpatialWeights(cds, reduction_method = reduction_method)
  lg <- calculateLocalG(cds = cds, lw = spatial_weights$lw,
                        wc = spatial_weights$wc, column_name = compare_col, var = compare)
  lg.df <- mergeLocalG(lg, df, compare_col) %>% getPval(method, tail)
  lg.df
}

### Functions for differential cell abundance testing across individuals. 

#' @title Convert a gene expression cell data set into a cell count cell data set.
#' @description Convert a gene expression cell data set into a cell count cell data set.
#'
#' @param cds (cell_data_set): Monocle3 cell data set object with gene expression information. 
#' @param sample_group (character): Column of cell metadata to group individual samples (e.g. replicates).
#' @param cell_group (character): Column of cell metadata to group counts.
#' @param sample_metadata (tibble): Additional metadata to add for sample data.
#' @param cell_metadata (tibble): Additional metadata to add for cell data.
#' @param lower_threshold (numeric): Minimum number of cells allowed for a group (e.g. cell type).
#' @param upper_threshold (numeric): Maximum number of cells allowed for a group (e.g. cell type).
#' @param covariate_cols (character): columns to include as covariates.  
#'
#' @examples
#' #'\dontrun{
#' cell_cds = make_cell_count_cds(cds, 
#' sample_group = "embryo", 
#' cell_group = "cell_type")
#'}
#' @return (cell_data_set): Cell_data_set containing a cell type compositions instead of gene expression data.
#'
#' @importFrom dplyr %>%
#' @export
make_cell_count_cds <- function(cds, 
                                sample_group,
                                cell_group, 
                                sample_metadata = NULL,
                                cell_metadata = NULL,
                                lower_threshold = NULL,
                                upper_threshold = NULL,
                                covariate_cols=NULL) {
  
  # check if anything contains NAs in it
  # if so drop them
  num_sample_group_NAs = sum(is.na(colData(cds)[[sample_group]]))
  if (num_sample_group_NAs != 0) {
    cds = cds[, !is.na(colData(cds)[[sample_group]])]
  }
  num_cell_group_NAs = sum(is.na(colData(cds)[[cell_group]]))
  if (num_cell_group_NAs != 0) {
    cds = cds[, !is.na(colData(cds)[[cell_group]])]
  }
  
  colData(cds)$sample <- NULL
  coldata_df <- colData(cds) %>% as.data.frame()
  
  coldata_df <- coldata_df %>% dplyr::rename("sample" = sample_group,
                                             "cell_group" = as.character(cell_group))
  
  coldata_df <- coldata_df %>% 
    group_by(sample, cell_group) %>% 
    mutate(group_id = cur_group_id()) %>% 
    ungroup()
  
  colData(cds)$group_id <- coldata_df$group_id
  
  cds_summary <- coldata_df %>%
    dplyr::group_by(sample, cell_group) %>%
    dplyr::summarize(cells = dplyr::n())
  
  cds_covariates_df <- coldata_df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarize(across(where(is.numeric), mean),
                     across(where(is.factor), function(x) { tail(names(sort(table(x))), 1) }),
                     across(where(is.character), function(x) { tail(names(sort(table(x, useNA="ifany"))), 1) } ))
  
  if (is.null(covariate_cols) == FALSE) {
    cds_covariates_df = cds_covariates_df %>% dplyr::select(sample, group_id, covariate_cols)
  }
  
  if (is.null(sample_metadata) == FALSE){
    cds_covariates_df = left_join(cds_covariates_df, sample_metadata, by=c("sample"="sample"))
    
  }
  
  cds_covariates_df <- cds_covariates_df %>% as.data.frame(cds_covariates_df)
  row.names(cds_covariates_df) <- cds_covariates_df %>% dplyr::pull(sample)
  
  cell_counts_wide <- tidyr::pivot_wider(cds_summary, names_from = sample, 
                                         values_from = cells, 
                                         values_fill = c(0))
  cell_states <-  as.character(cell_counts_wide %>% dplyr::pull("cell_group"))
  cell_counts_wide <-  as.matrix(cell_counts_wide[,2:ncol(cell_counts_wide)])
  
  row.names(cell_counts_wide) = cell_states
  
  # filter out cell groups based on counts
  if (is.null(lower_threshold) == FALSE) {
    cell_counts_wide = cell_counts_wide[Matrix::rowSums(cell_counts_wide) >= lower_threshold,]
  }
  if (is.null(upper_threshold) == FALSE) {
    cell_counts_wide = cell_counts_wide[Matrix::rowSums(cell_counts_wide) <= upper_threshold,]
  }
  
  cds_covariates_df = cds_covariates_df[colnames(cell_counts_wide),]
  
  cell_count_cds = monocle3::new_cell_data_set(cell_counts_wide,
                                               cell_metadata=cds_covariates_df,
                                               gene_metadata=cell_metadata)
  
  return(cell_count_cds)
}



#' @title Convert a gene expression cell data set into a cell count cell data set.
#' @description Convert a gene expression cell data set into a cell count cell data set.
#'
#' @param ccs (cell_data_set): Cell count cell_data_set. 
#' @param comp_col (character): Column of cell metadata for comparing cell type abundances.
#' @param model_formula (character): Model variable for beta binomial regression testing.
#' @param sample_metadata (tibble): Additional metadata to add for sample data.
#' @param cell_metadata (tibble): Additional metadata to add for cell data.
#' @param lower_threshold (numeric): Minimum number of cells allowed for a group (e.g. cell type).
#' @param upper_threshold (numeric): Maximum number of cells allowed for a group (e.g. cell type).
#' @param covariate_cols (character): columns to include as covariates.  
#'
#' @examples
#'\dontrun{
#'bb_res = bb_compare_abundance(cell_cds, 
#'                              ctrl_ids = c("control_genotype", "test_genotype"),
#'                              comp_col = "gene_target",
#'                              model_formula = "~ gene_target")
#'}
#'
#' @return (cell_data_set): Cell_data_set containing a cell type compositions instead of gene expression data.
#'
#' @importFrom dplyr %>%
#' @export
bb_compare_abundance <- function(ccs,
                                 comp_col,
                                 model_formula, 
                                 ctrl_ids,
                                 nuisance_cols = NULL,
                                 nuisance_formula = NULL,
                                 ...) {
  
  # make ccs
  colData(ccs)[[comp_col]] = as.data.frame(colData(ccs)) %>%
    mutate(comp_group = ifelse(!!sym(comp_col) %in% ctrl_ids, "ctrl", !!sym(comp_col))) %>%
    pull(comp_group)
  
  ccs_coldata = colData(ccs) %>% as.data.frame() %>% dplyr::select(sample, !!sym(comp_col), 
                                                                   Size_Factor, nuisance_cols) %>% 
    ungroup()
  
  count_df = SingleCellExperiment::counts(ccs) %>%
    as.matrix() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_group") %>%
    tidyr::pivot_longer(cols = !cell_group, names_to = "sample", values_to = "cells") %>%
    dplyr::group_by(sample) %>%
    left_join(ccs_coldata, by="sample") %>%
    dplyr::mutate(cells = round(cells/Size_Factor)) %>%
    dplyr::mutate(total_cells = sum(cells))
  
  fc_df = count_df %>%
    dplyr::group_by_at(vars(comp_col, cell_group, nuisance_cols)) %>%
    dplyr::summarize(cell_mean = mean(cells)) %>%
    tidyr::pivot_wider(names_from = comp_col, values_from = "cell_mean") %>%
    tidyr::pivot_longer(-c(cell_group, ctrl, nuisance_cols), names_to = comp_col, values_to = "cells") %>%
    dplyr::mutate(abund_log2fc = log2((cells + 1)/(ctrl+1))) %>%
    dplyr::arrange(!!sym(comp_col))
  
  
  # to fix new col name
  count_df[[comp_col]] = forcats::fct_relevel(count_df[[comp_col]], "ctrl")
  
  cell.groups = unique(count_df$cell_group)
  
  # a help function to tidy the vgam model output - used in compare_abundance function
  tidy.vglm = function(x, conf.int=FALSE, conf.level=0.95) {
    co <- as.data.frame(coef(summary(x)))
    names(co) <- c("estimate","std.error","statistic","p.value")
    if (conf.int) {
      qq <- qnorm((1+conf.level)/2)
      co <- transform(co,
                      conf.low=estimate-qq*std.error,
                      conf.high=estimate+qq*std.error)
    }
    co <- data.frame(term=rownames(co),co)
    rownames(co) <- NULL
    return(co)
  }
  
  fit_beta_binomial = function(cg, count_df, model_formula, ...) {
    type_df = count_df %>% filter(cell_group == cg)
    count_df = cbind(type_df$cells, type_df$total_cells - type_df$cells)
    fit =  VGAM::vglm(as.formula(model_formula), betabinomial, data = type_df, trace = TRUE, ...)
    fit_df = tidy.vglm(fit)
    
  }
  
  model_formula = stringr::str_replace_all(model_formula, "~", "")
  model_formula_str = paste0("count_df ~", model_formula)
  
  if (is.null(nuisance_formula) == FALSE) {
    nuisance_formula = stringr::str_replace_all(nuisance_formula, "~", "")
    model_formula_str = paste0(model_formula_str, "+", nuisance_formula)
  }
  
  bb_res = data.frame("cell_group" = cell.groups) %>% 
    dplyr::mutate("beta_binomial" = purrr::map(.f  = purrr::possibly(fit_beta_binomial, NA_character_),
                                               .x = cell_group,
                                               count_df = count_df,
                                               model_formula = model_formula_str)) %>%
    dplyr::filter(!is.na(beta_binomial)) %>%
    tidyr::unnest(c(beta_binomial)) %>%
    dplyr::arrange(desc(estimate)) %>%
    dplyr::filter(!grepl("Intercept", term))
  
  bb_res = left_join(bb_res,
                     fc_df %>% dplyr::select(cell_group, abund_log2fc, "ctrl_mean" = "ctrl"),
                     by = "cell_group") 
  
  return(bb_res)
  
}

