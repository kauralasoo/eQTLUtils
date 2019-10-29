#' Generate distance matrix of genotype samples using PC matrices
#'
#' @param reference_pca_df PC matrix of reference genotype data. Should have genotype_id, superpopulation_code and PCN columns
#' @param pca_df PC matrix of input genotype data. Should have genotype_id, superpopulation_code and PCN columns
#' @param n_pcs number of PCs to be taken into account. Default value is 3
#' @param method Method of distance measure. See stats::dist() method for options
#' @author Nurlan Kerimov
#' @export
generate_distance_matrix <- function(reference_pca_df, pca_df, n_pcs = 3, method = "euclidean"){
  # required fields assertion ====
  assertthat::assert_that(hasName(pca_df, "genotype_id"), msg = "Column genotype_id is missing in pca_df")
  assertthat::assert_that(hasName(pca_df, "superpopulation_code"), msg = "Column superpopulation_code is missing in pca_df")
  pcs <- c()
  for (i in 1:n_pcs) {
    pc_num <- paste0("PC",i)
    assertthat::assert_that(hasName(pca_df, pc_num), msg = paste0("Column ", pc_num, " is missing in pca_df"))
    pcs <- c(pcs, pc_num)
  }
  
  assertthat::assert_that(hasName(reference_pca_df, "genotype_id"), msg = "Column genotype_id is missing in reference_pca_df")
  assertthat::assert_that(hasName(reference_pca_df, "superpopulation_code"), msg = "Column superpopulation_code is missing in reference_pca_df")
  for (i in 1:n_pcs) {
    pc_num <- paste0("PC",i)
    assertthat::assert_that(hasName(reference_pca_df, pc_num), msg = paste0("Column ", pc_num, " is missing in reference_pca_df"))
  }
  # ====
  
  sample_size <- pca_df %>% nrow()
  reference_pca_df$genotype_id <- paste0("ref_", reference_pca_df$genotype_id)
  comb_pcs = rbind(pca_df, reference_pca_df)
  rownames(comb_pcs) <- comb_pcs$genotype_id
  distance_matrix <- dist(comb_pcs[pcs], method = method) %>% as.matrix(labels = TRUE)
  distance_matrix <- distance_matrix[-c(1:sample_size), 1:sample_size] 
  tibble_ind <- distance_matrix %>% as_tibble(rownames="genotype_id")
  tibble_ind <- reference_pca_df[,c("genotype_id","superpopulation_code")] %>% left_join(tibble_ind)
  
  mean_ind_diff_summary <- tibble_ind[,-1] %>% group_by(superpopulation_code) %>% summarize_all("mean")
  return(mean_ind_diff_summary)
}

#' Assigns populations using distance matrix from generate_distance_matrix() function
#'
#' @param distance_matrix distance matrix dataframe created by generate_distance_matrix() function
#' @param abs_threshold sample's minimum distance to any reference population should be less than this threshold to be assigned to population. Default value is 0.02
#' @param rel_threshold sample's minimum distance to any reference population should be rel_threshold times smaller than second minimum distanceto be assigned to population. Default value is 1.7
#' @author Nurlan Kerimov
#' @export
assign_populations <- function(distance_matrix, abs_threshold = 0.02, rel_threshold = 1.7){
  assigned <- distance_matrix[,-1] %>% 
    rbind(if_else(apply(distance_matrix[,-1], 2, min)<abs_threshold, 
                  as.character(distance_matrix$superpopulation_code[apply(distance_matrix[,-1], 2, which.min)]), 
                  "Admixed")) %>% as.data.frame()
  rownames(assigned) <- c(as.character(distance_matrix$superpopulation_code), "pop_assign_abs_thresh")
  assigned <- assigned %>% t() %>% as.data.frame() 
  assigned <- cbind(genotype_id=rownames(assigned), assigned)
  
  # find if there are admixed samples
  mins <- apply(distance_matrix[,-1], 2, sort) %>% t() %>% as.data.frame() 
  mins <- mins %>% mutate(isAdmixed = mins[,1]*rel_threshold > mins[,2])
  # set assigned_population as minimum distance population first
  assigned_rel <- distance_matrix[,-1] %>% rbind(as.character(distance_matrix$superpopulation_code[apply(distance_matrix[,-1], 2, which.min)])) %>% as.data.frame()
  rownames(assigned_rel) <- c(as.character(distance_matrix$superpopulation_code), "pop_assign_rel_thresh")
  assigned_rel <- assigned_rel %>% t() %>% as.data.frame()
  
  # if there are admixed samples according to proportional threshold replace assigned_rel population value
  if (sum(mins$isAdmixed)>0) {
    levels(assigned_rel$pop_assign_rel_thresh) <- c(levels(assigned_rel$pop_assign_rel_thresh), "Admixed")
    assigned_rel$pop_assign_rel_thresh[mins$isAdmixed] <- "Admixed"
  }
  assigned_rel <- cbind(genotype_id=rownames(assigned_rel), assigned_rel)
  
  return(left_join(assigned, assigned_rel[,c("genotype_id", "pop_assign_rel_thresh")]))
}

#' Reassigns already assigned populations according to new thresholds
#'
#' @param distance_matrix distance matrix dataframe created by assign_populations() function
#' @param abs_threshold sample's minimum distance to any reference population should be less than this threshold to be assigned to population
#' @param rel_threshold sample's minimum distance to any reference population should be rel_threshold times smaller than second minimum distanceto be assigned to population
#' @author Nurlan Kerimov
#' @export
reassign_populations <- function(distance_matrix, abs_threshold, rel_threshold){
  distance_matrix <- distance_matrix[c(2:5)] %>% t() %>% as.data.frame() %>% setNames(distance_matrix$genotype_id)
  distance_matrix <- cbind(superpopulation_code = rownames(distance_matrix), distance_matrix)
  assigned_new_pop <- assign_populations(distance_matrix, abs_threshold = 0.01, rel_threshold = 3)
  return(assigned_new_pop)
}