estimateMedianTPM <- function(se, subset_by = "qtl_group", assay_name = "counts"){
  
  #Extract colData from SE
  col_data = SummarizedExperiment::colData(se)
  assertthat::assert_that(assertthat::has_name(col_data, "study"))
  assertthat::assert_that(assertthat::has_name(col_data, subset_by))
  
  #Calculate mean and median expression per QTL group
  tpm_se = normaliseSE_tpm(se, assay_name = "counts")
  
  #Extract study id
  study_id = col_data$study[1]
  
  #Split into groups
  groups = unique(col_data[,subset_by])
  group_list = setNames(groups, groups)
  group_se_list = purrr::map(qtl_group_list, ~subsetSEByColumnValue(tpm_se, subset_by, .))
  
  #Estimate median TPM per group
  median_list = purrr::map(group_se_list, ~SummarizedExperiment::assays(.)[["tpms"]] %>% apply(.,1,median))
  
  #Make into a single data frame
  median_df = purrr::map_df(median_list, ~dplyr::data_frame(phenotype_id = names(.), median_tpm = .), .id = subset_by)
  
  #reorder results
  results_df = dplyr::mutate(median_df, study = study_id) %>%
    dplyr::select(phenotype_id, study, everything()) %>%
    dplyr::mutate(median_tpm = round(median_tpm,3))
  return(results_df)
}