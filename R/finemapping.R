calculatePairwiseCSMetrics <- function(cs_list){
  
  #Find all combinations
  combinations = combn(length(cs_list),2, simplify = F)
  
  #Apply to pairs
  cs_pairs = purrr::map(combinations, ~cs_list[.])
  pairwise_similarities = purrr::map_df(cs_pairs, ~calculateSinglePairMetrics(.))
  return(pairwise_similarities)
}

calculateSinglePairMetrics <- function(single_pair){
  metrics = dplyr::tibble(cs1 = names(single_pair)[1], 
                          cs2 = names(single_pair)[2], 
                          cs1_size = length(single_pair[[1]]),
                          cs2_size = length(single_pair[[2]]),
                          intersect = length(intersect(single_pair[[1]], single_pair[[2]])),
                          union = length(union(single_pair[[1]], single_pair[[2]]))
  ) %>%
    dplyr::mutate(jaccard = intersect/union,
                  max_subset = intersect/min(cs1_size, cs2_size))
  return(metrics)
}


importSusieCredibleSets <- function(cs_folder_path){
  
  #Define file paths
  file_names = list.files(cs_folder_path)
  full_paths = file.path(cs_folder_path, file_names)
  path_list = setNames(full_paths, file_names)
  
  #Import files
  cs_df = purrr::map_df(path_list, ~readr::read_tsv(.), .id = "file_name") %>%
    tidyr::separate(file_name, c("study","qtl_group", "quant_method", "chunk", "suffix"), sep = "\\.") %>%
    dplyr::select(-chunk, -suffix)
  return(cs_df)
}
