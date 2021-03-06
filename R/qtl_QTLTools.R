#' Save a list of matrices into a suitable format for QTLTools
#'
#' Works with expression and covariates matrices.
#'
#' @param data_list list of matrices
#' @param output_dir relative path to the output dir
#' @param file_suffix suffix added to each file after their name in the list.
#' @param file_suffix suffix added to each file after their name in the list.
#' @return None
#' @author Kaur Alasoo
#' @export
saveQTLToolsMatrices <- function(data_list, output_dir, file_suffix = "bed", file_prefix = "", col_names = TRUE){

  #Check if the output dir exists and if not then create one
  if(!file.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }

  #Save each matrix as a separate  txt file
  for (sn in names(data_list)){
    file_name = ifelse(file_prefix == "", paste(sn, file_suffix, sep = "."), paste(file_prefix, sn, file_suffix, sep = "."))
    file_path = file.path(output_dir,file_name)
    print(file_path)
    write.table(data_list[[sn]], file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = col_names)
  }
}

#' Import QTLtools output table from permutation run into R.
#'
#'
#' @param file_path Path to the QTLtools output file.
#' @return data_frame containing gene_ids, snp ids and p-values.
#' @author Kaur Alasoo
#' @export
importQTLtoolsTable <- function(file_path){
  col_names = c("group_id", "pheno_chr", "pheno_start", "pheno_end", "strand", "phenotype_id", "group_size", "n_cis_snps",
                "distance", "snp_id", "snp_chr", "snp_start", "snp_end", "df", "dummy", "beta1",
                "beta2", "p_nominal","slope","p_perm","p_beta")
  col_types = "cciicciiicciiiddddddd"
  table = readr::read_delim(file_path, col_names = col_names, delim = " ", col_types = col_types) %>%
    dplyr::filter(!is.na(p_beta)) %>%
    dplyr::mutate(p_bonferroni = p_nominal*group_size*n_cis_snps) %>%
    dplyr::mutate(p_bonferroni = pmin(p_bonferroni,1)) %>%
    dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>%
    dplyr::arrange(p_fdr)
  return(table)
}

#' Import QTLtools output table from permutation run v2 into R.
#'
#'
#' @param file_path Path to the QTLtools permutation output file.
#' @return data_frame containing gene_ids, snp ids and p-values.
#' @author Nurlan Kerimov
#' @export
importQTLtoolsPermRes<- function(file_path){
  col_types = "cciiciddddd"
  table = readr::read_tsv(file_path, col_types = col_types) %>%
    dplyr::filter(!is.na(p_beta)) %>%
    dplyr::mutate(p_bonferroni = pvalue*n_variants*n_traits) %>%
    dplyr::mutate(p_bonferroni = pmin(p_bonferroni,1)) %>%
    dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>%
    dplyr::arrange(p_fdr)
  return(table)
}

importQTLtoolsNominalTable <- function(file_path){
  col_names = c("phenotype_id", "pheno_chr", "pheno_start", "pheno_end", "strand", "n_cis_snps", "distance", "snp_id", "snp_chr", "snp_start", "snp_end", "p_nominal", "beta", "is_lead")
  col_types = "cciiciicciiddi"
  table = readr::read_tsv(file_path, col_names = col_names, col_types = col_types) %>%
    dplyr::group_by(phenotype_id) %>%
    dplyr::mutate(p_bonferroni = p.adjust(p_nominal, n = n_cis_snps)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, method = "fdr"))
  return(table)
}

#' Fetch particular genes from tabix indexed QTLtools output file.
#'
#' @param gene_ranges GRanges object with coordinates of the cis regions around genes.
#' @param tabix_file Tabix-indexed fastqtl output file.
#'
#' @return List of data frames containing QTLtools results for each gene.
#' @export
qtltoolsTabixFetchPhenotypes <- function(phenotype_ranges, tabix_file){

  #Assertions about input
  assertthat::assert_that(class(phenotype_ranges) == "GRanges")
  assertthat::assert_that(assertthat::has_name(GenomicRanges::elementMetadata(phenotype_ranges), "phenotype_id"))

  #Set column names for rasqual
  qtltools_columns = c("phenotype_id","pheno_chr","pheno_start", "pheno_end",
                      "strand","n_snps", "distance", "snp_id", "snp_chr",
                      "snp_start", "snp_end", "p_nominal","beta", "is_lead")
  qtltools_coltypes = "cciiciicciiddi"

  result = list()
  for (i in seq_along(phenotype_ranges)){
    selected_phenotype_id = phenotype_ranges[i]$phenotype_id
    print(i)
    tabix_table = scanTabixDataFrame(tabix_file, phenotype_ranges[i],
                                     col_names = qtltools_columns, col_types = qtltools_coltypes)[[1]]

    #Filter by phenotype id if the tabix table is not null
    if(!is.null(tabix_table)){
      tabix_table = dplyr::filter(tabix_table, phenotype_id == selected_phenotype_id)
    }

    #Add additional columns
    result[[i]] = tabix_table
  }
  return(result)
}

# Fetch phenotype_id and snp_id pairs from the tabix file. Extract variant coordinates
# directly from the snp_id column of the selected_pairs data frame.
qtltoolsTabixFetchPhenotypeVariantPairs <- function(selected_pairs, tabix_file){

  #Extract variant coordinates from variant id
  var_coords = tidyr::separate(selected_pairs, snp_id, into = c("chr", "pos", "ref", "alt"), remove = F) %>%
    dplyr::mutate(chr = stringr::str_replace(chr, "chr", "")) %>%
    dplyr::mutate(pos = as.integer(pos)) %>%
    dplyr::select(-ref, -alt)

  #Make a Granges object
  ranges = GenomicRanges::GRanges(seqnames = var_coords$chr,
                                  ranges = IRanges(start = var_coords$pos, end = var_coords$pos),
                                  strand = "*")
  elementMetadata(ranges) = dplyr::select(var_coords, -chr, -pos)

  #Fetch data from tabix
  tabix_data = eQTLUtils::qtltoolsTabixFetchPhenotypes(
    tabix_file = tabix_file,
    phenotype_ranges = ranges) %>%
    purrr::map_df(identity) %>%
    dplyr::semi_join(var_coords, by = c("phenotype_id", "snp_id"))
  return(tabix_data)
}


#' Post-process QTLTools mbv results to find the best matching individual for each sample
#'
#' @param mbv_df Data frame with MBV results for one sequencing sample.
#'
#' @return Data frame with one row identifying the best matching invidivual for this sample
#' @export
mbvFindBestMatch <- function(mbv_df){
  res = dplyr::transmute(mbv_df, mbv_genotype_id = SampleID,
                         het_consistent_frac = n_het_consistent/n_het_covered,
                         hom_consistent_frac = n_hom_consistent/n_hom_covered)

  #Identify best het
  best_het = dplyr::arrange(res, -het_consistent_frac) %>% dplyr::filter(dplyr::row_number() == 1)
  other_het = dplyr::arrange(res, -het_consistent_frac) %>% dplyr::filter(dplyr::row_number() > 1)
  best_row = dplyr::mutate(best_het, het_min_dist = min(best_het$het_consistent_frac - other_het$het_consistent_frac),
                           hom_min_dist = min(best_het$hom_consistent_frac - other_het$hom_consistent_frac),
                           distance = sqrt(het_min_dist^2 + hom_min_dist^2))

  #Compare against best hom
  best_hom = dplyr::arrange(res, -hom_consistent_frac) %>% dplyr::filter(dplyr::row_number() == 1)
  if(best_row$mbv_genotype_id != best_hom$mbv_genotype_id){
    best_row = dplyr::mutate(best_row, het_consistent_frac = as.numeric(NA), hom_consistent_frac = as.numeric(NA),
                             het_min_dist = as.numeric(NA), hom_min_dist = as.numeric(NA))
  }
  return(best_row)
}


#' Import MBV output files for all samples in a studt
#'
#' @param mbv_dir directory where mbv output files are located
#'
#' @return List of the matrixes with mbv output values per sample
#' @export
mbvImportData <- function(mbv_dir, suffix = ".mbv_output.txt"){

  #List all files
  mbv_files = list.files(mbv_dir, full.names = T)
  mbv_files = mbv_files[grep(suffix, mbv_files)]

  #Make sample names
  sample_names = stringr::str_replace_all(basename(mbv_files), suffix, "")
  sample_list = setNames(mbv_files, sample_names)

  #Import mbv files
  mbv_results = purrr::map(sample_list, ~readr::read_delim(., delim = " ", col_types = "ciiiiiiiiii"))

  return(mbv_results)
}

#' Convert SummarizedExperiment object into a bed file suitable for QTLTools
#'
#' @param se SummarizedExperiment object
#' @param assay_name Assay of the SummarizedExperiment object used for QTL mapping
#'
#' @return Assay converted into bed format suitable for QTLTools
#' @export
convertSEtoQTLtools <- function(se, assay_name = "cqn"){

  #Extract rowData from the SE
  phenotype_data = rowData(se) %>%
    as.data.frame() %>%
    dplyr::as_tibble()

  #Make sure that all required columns are present
  assertthat::assert_that(assertthat::has_name(phenotype_data, "chromosome"))
  assertthat::assert_that(assertthat::has_name(phenotype_data, "phenotype_pos"))
  assertthat::assert_that(assertthat::has_name(phenotype_data, "phenotype_id"))
  assertthat::assert_that(assertthat::has_name(phenotype_data, "group_id"))
  assertthat::assert_that(!is.null(se$genotype_id))

  #Make genePos table for QTLTools
  pheno_data = dplyr::arrange(phenotype_data, chromosome, phenotype_pos) %>%
    dplyr::transmute(chromosome, left = phenotype_pos, right = phenotype_pos, phenotype_id, group_id, strand) %>%
    dplyr::rename_("#chr" = "chromosome") %>%
    dplyr::mutate(strand = ifelse(strand == 1, "+", "-"))

  #Exptract phenotype and rename columns according to genotype id
  assay = assays(se)[[assay_name]]
  colnames(assay) = se$genotype_id
  assay = round(assay, 3) #Round to three digits

  #Make QTLtools phenotype table
  res = dplyr::mutate(as.data.frame(assay), phenotype_id = rownames(assay)) %>%
    dplyr::select(phenotype_id, dplyr::everything()) %>%
    dplyr::left_join(pheno_data, ., by = "phenotype_id") %>%
    dplyr::arrange()

  return(res)
}


importQTLtoolsPCA <- function(pca_path){
  naive_pca = readr::read_delim(pca_path, delim = " ") %>%
    dplyr::select(-SampleID)
  sample_ids = colnames(naive_pca)
  pca_df = t(naive_pca) %>%
    as.data.frame() %>%
    dplyr::as_tibble()
  colnames(pca_df) = paste0("PC", 1:ncol(pca_df))
  pca_df$sample_id = sample_ids
  pca_df = dplyr::select(pca_df, sample_id, everything())
  return(pca_df)
}

#' Split study SummarizedExperiment object into invidiual input files for QTLTools.
#'
#' Uses the qtl_group metadata column (required) to split a single SummarizedExperiment
#' object into multiple independent input files for QTLTools. Write the files directly to disk.
#'
#' @param se SummarizedExperiment object.
#' @param assay_name Name of the assay in the SummarizedExperiment object used for QTL mapping.
#' @param out_dir Path to the output directory where the QTLTools input files will be written.
#' @param extra_qtl_group Used to split the datasets using some additional metadata column in additon
#' to the qtl_group column. This can for example be used to perform sex-specific or population-specific
#' eQTL analysis.
#'
#' @return None
#' @export
studySEtoQTLTools <- function(se, assay_name, out_dir, extra_qtl_group = NULL){

  #Make assertions
  assertthat::assert_that(assertthat::has_name(SummarizedExperiment::colData(se), "qtl_group"))
  assertthat::assert_that(assertthat::has_name(SummarizedExperiment::assays(se), assay_name))

  if(is.null(extra_qtl_group)){
    #Split the SE into list based on qtl_group
    qtl_groups = unique(se$qtl_group)
    group_list = setNames(as.list(qtl_groups), qtl_groups)
    group_se_list = purrr::map(group_list, ~subsetSEByColumnValue(se, "qtl_group", .))

   } else{

    #If extra_qtl_group is specified, then split the dataset using both the qtl_group values
    #as well as the values of the extra_qtl_group column.
    assertthat::assert_that(assertthat::has_name(SummarizedExperiment::colData(se), extra_qtl_group))
    col_data = SummarizedExperiment::colData(se)
    se$new_qtl_group = paste(col_data[,"qtl_group"], col_data[,extra_qtl_group],sep = "_")
    qtl_groups = unique(se$new_qtl_group)
    group_list = setNames(as.list(qtl_groups), qtl_groups)
    group_se_list = purrr::map(group_list, ~subsetSEByColumnValue(se, "new_qtl_group", .))
  }

  #Convert SE onbjects to QTLtools
  qtltools_list = purrr::map(group_se_list, ~convertSEtoQTLtools(., assay_name = assay_name))
  saveQTLToolsMatrices(qtltools_list, output_dir = out_dir, file_suffix = "bed")

  #Extract sample names
  sample_names = purrr::map(qtltools_list, ~colnames(.)[-(1:6)])
  saveQTLToolsMatrices(sample_names, output_dir = out_dir, file_suffix = "sample_names.txt", col_names = FALSE)
  
  count_matrices <- purrr::map(group_se_list, ~SummarizedExperiment::cbind(phenotype_id = rownames(assays(.)[[assay_name]]), assays(.)[[assay_name]]))
  saveQTLToolsMatrices(count_matrices, output_dir = output_dir, file_suffix = "tsv")
}

#' Split study SummarizedExperiment object into invidiual input files for QTLTools.
#'
#' Uses the qtl_group metadata column (required) to split a single SummarizedExperiment
#' object into multiple independent input files for QTLTools. Write the files directly to disk.
#'
#' @param se SummarizedExperiment object.
#' @param assay_name Name of the assay in the SummarizedExperiment object used for QTL mapping.
#' @param out_dir Path to the output directory where the QTLTools input files will be written.
#' @param study_name Custom study_name character string to safe files with prefix
#'
#' @return None
#' @export
studySEtoCountMatrices <- function(se, assay_name, out_dir, study_name = NULL, quantile_tpms = NULL, tpm_thres = 0.1){
  #Make assertions
  assertthat::assert_that(assertthat::has_name(SummarizedExperiment::colData(se), "qtl_group"))
  assertthat::assert_that(assertthat::has_name(SummarizedExperiment::assays(se), assay_name))
  
  if (base::is.null(study_name)) {
    study_name = unique(colData(se)["study"])[[1]]
  }
  
  #Split the SE into list based on qtl_group
  qtl_groups = unique(se$qtl_group)
  group_list = setNames(as.list(qtl_groups), qtl_groups)
  group_se_list = purrr::map(group_list, ~subsetSEByColumnValue(se, "qtl_group", .))
  if (!is.null(quantile_tpms)) {
    group_se_list <- purrr::map(group_se_list, ~filterTPMQuantile(., quantile_tpms = quantile_tpms, tpm_thres = tpm_thres))
  }
  count_matrices <- purrr::map(group_se_list, ~SummarizedExperiment::cbind(phenotype_id = rownames(assays(.)[[assay_name]]), assays(.)[[assay_name]]))
  saveQTLToolsMatrices(count_matrices, output_dir = file.path(out_dir, "qtl_group_split_norm") , file_suffix = "tsv", file_prefix = study_name)
}

filterTPMQuantile <- function(subset_se, quantile_tpms, tpm_thres = 1){
  #Check that required columns exist
  assertthat::assert_that(assertthat::has_name(quantile_tpms, "qtl_group"))
  assertthat::assert_that(assertthat::has_name(quantile_tpms, "median_tpm"))
  assertthat::assert_that(assertthat::has_name(quantile_tpms, "phenotype_id"))
  
  #Find expressed genes
  sample_meta_qtlgroup <- SummarizedExperiment::colData(subset_se) %>% SummarizedExperiment::as.data.frame()
  phenotype_data <- SummarizedExperiment::rowData(subset_se) %>% SummarizedExperiment::as.data.frame()
  assertthat::assert_that(assertthat::assert_that(sample_meta_qtlgroup$qtl_group %>% unique() %>% length()==1), 
                          msg = "There are more than 1 qtl_groups in qtlgroup subset")
  selected_qtl_group = sample_meta_qtlgroup$qtl_group %>% unique()
  message("Filter SE by quntile TPMs for QTL group: ", selected_qtl_group)
  
  not_expressed_genes = dplyr::filter(quantile_tpms, qtl_group == selected_qtl_group, median_tpm < tpm_thres)
  
  #Find expressed phenotyes
  expressed_phenotypes = setdiff(phenotype_data$gene_id, not_expressed_genes$phenotype_id)
  message(paste0("Number of expressed genes included in the analysis: ", length(expressed_phenotypes)))
  expressed_phenotype_metadata = dplyr::filter(phenotype_data, gene_id %in% expressed_phenotypes)
  
  return (subset_se[rowData(subset_se)$phenotype_id %in% expressed_phenotype_metadata$phenotype_id,])
}
