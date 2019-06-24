message(" ## Loading libraries: optparse")
suppressPackageStartupMessages(library("optparse"))

#Parse command-line options
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  make_option(c("-c", "--count_matrix"), type="character", default=NULL,
              help="Counts matrix file path. Tab separated file", metavar = "type"),
  make_option(c("-s", "--sample_meta"), type="character", default=NULL,
              help="Sample metadata file. Tab separated file", metavar = "type"),
  make_option(c("-p", "--phenotype_meta"), type="character", default=NULL,
              help="Phenotype metadata file. Tab separated file", metavar = "type"),
  make_option(c("-q", "--quant_method"), type="character", default="ge",
              help="Quantification method. Possible values: ge, tx, ex and sp. [default \"%default\"]", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default="./RNA_QC_RESULTS/",
              help="Path to the output directory. [default \"%default\"]", metavar = "type"),
  make_option(c("-g", "--generate_plots"), type="logical", default=TRUE,
              help="Flag to generate the plots. [default \"%default\"]", metavar = "bool"),
  make_option(c("--build_html"), type="logical", default=FALSE,
              help="Flag to build plotly html plots [default \"%default\"]", metavar = "bool"),
  make_option(c("-m", "--mbvdir"), type="character", default=NULL,
              help="Path to the output directory. [default \"%default\"]", metavar = "type"),
  make_option(c("-n", "--name_of_study"), type="character", default=NULL,
              help="Name of the study. Optional", metavar = "type")
)

message(" ## Parsing options")
opt <- parse_args(OptionParser(option_list=option_list))

message(" ## Loading libraries: devtools, dplyr, SummarizedExperiment, cqn, data.table")
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("cqn"))
suppressPackageStartupMessages(library("data.table"))

#Debugging
if (FALSE) {
  opt = list()
  opt$c="data/featureCounts_matrices/Alasoo_test_data/Alasoo_merged_gene_counts.txt"
  opt$s="data/sample_metadata/Alasoo_2018.tsv"
  opt$p="data/annotations/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz"
  opt$m="mbv/"
}

count_matrix_path = opt$c
sample_meta_path = opt$s
phenotype_meta_path = opt$p
output_dir = opt$o
generate_plots = opt$g
build_html = opt$build_html
mbv_files_dir = opt$m
quant_method = opt$q
study_name = opt$n

message("######### Options: ######### ")
message("######### Working Directory  : ", getwd())
message("######### quant_method       : ", quant_method)
message("######### count_matrix_path  : ", count_matrix_path)
message("######### sample_meta_path   : ", sample_meta_path)
message("######### phenotype_meta_path: ", phenotype_meta_path)
message("######### output_dir         : ", output_dir)
message("######### generate_plots     : ", generate_plots)
message("######### build_html         : ", build_html)
message("######### mbv_files_dir      : ", mbv_files_dir)
message("######### opt_study_name     : ", study_name)

if (!dir.exists(paste0(output_dir, "/normalised/"))){
  dir.create(paste0(output_dir, "/normalised/"), recursive = TRUE)
}

if (build_html) { 
  message(" ## Loading libraries: plotly")
  suppressPackageStartupMessages(library("plotly")) 
}

# Read the inputs
message("## Reading sample metadata ##")
sample_metadata <- utils::read.csv(sample_meta_path, sep = '\t')

message("## Reading featureCounts transcript metadata ##")
phenotype_meta = readr::read_delim(phenotype_meta_path, delim = "\t", col_types = "ccccciiicciidi")

message("## Reading featureCounts matrix ##")
data_fc <- utils::read.csv(count_matrix_path, sep = '\t')

message("## Make Summarized Experiment ##")
se <- eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = data_fc, row_data = phenotype_meta, col_data = sample_metadata)

if (!dir.exists(paste0(output_dir, "/rds/"))){
  dir.create(paste0(output_dir, "/rds/"), recursive = TRUE)
}
if (!dir.exists(paste0(output_dir, "/tsv/"))){
  dir.create(paste0(output_dir, "/tsv/"), recursive = TRUE)
}

#add assertion checks for needed columns
if (is.null(study_name)) { 
  assertthat::has_name(sample_metadata, "study" )
  study_name <- sample_metadata$study[1] 
}

message("## Perform PCA calculation ##")
pca_res <- eQTLUtils::plotPCAAnalysis(study_data_se = se, export_output = generate_plots, html_output = build_html, output_dir = output_dir)
saveRDS(pca_res, paste0(output_dir, paste0("/rds/", study_name ,"_pca_res.rds")))
readr::write_tsv(pca_res$pca_matrix, paste0(output_dir, paste0("/tsv/", study_name ,"_pca_matrix.tsv")))

message("## Perform MDS calculation ##")
mds_res <- eQTLUtils::plotMDSAnalysis(study_data_se = se, export_output = generate_plots, html_output = build_html, output_dir = output_dir)
saveRDS(mds_res, paste0(output_dir, paste0("/rds/", study_name ,"_mds_res.rds")))
readr::write_tsv(mds_res, paste0(output_dir, paste0("/tsv/", study_name ,"_mds_matrix.tsv")))

message("## Perform MDS calculation ##")
sex_spec_gene_exp <- eQTLUtils::plotSexQC(study_data = se, export_output = generate_plots, html_output = build_html, output_dir = output_dir)
saveRDS(sex_spec_gene_exp, paste0(output_dir, paste0("/rds/", study_name ,"_sex_spec_gene_exp_res.rds")))
readr::write_tsv(sex_spec_gene_exp, paste0(output_dir, paste0("/tsv/", study_name ,"_sex_spec_gene_exp_matrix.tsv")))

if (!is.null(mbv_files_dir)) {
  message("## Perform MBV Analysis ##")
  mbv_results = eQTLUtils::mbvImportData(mbv_dir = mbv_files_dir, suffix = ".mbv_output.txt")
  best_matches = purrr::map_df(mbv_results, eQTLUtils::mbvFindBestMatch, .id = "sample_id") %>% dplyr::arrange(distance)
  best_matches <- left_join(sample_metadata[,c("sample_id","genotype_id")], best_matches)
  best_matches$is_correct_match <- best_matches$mbv_genotype_id == best_matches$genotype_id
  readr::write_tsv(best_matches, paste0(output_dir, paste0("/tsv/", study_name ,"_MBV_best_matches_matrix.tsv"))) 
  
  eQTLUtils::plot_mbv_results(mbv_files_path = mbv_files_dir, output_path = paste0(output_dir, "/MBV/"))
}
message("## RNA Quality Control is completed! ##")


















