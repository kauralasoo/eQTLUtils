makeFeatureCountsSummarizedExperiment <- function(featureCounts, transcript_metadata, sample_metadata){

  #Specify required gene metadata columns
  required_gene_meta_columns = c("phenotype_id","quant_id","group_id","gene_id","chromosome","gene_start",
                                 "gene_end","strand","gene_name","gene_type","phenotype_gc_content","gene_version","phenotype_pos")

  #Extract gene metadata
  gene_data = dplyr::select(transcript_metadata, gene_id, chromosome, gene_start, gene_end, strand,
                            gene_name, gene_type, phenotype_gc_content, gene_version) %>%
    dplyr::distinct() %>%
    dplyr::mutate(phenotype_id = gene_id, group_id = gene_id, quant_id = gene_id) %>%
    dplyr::mutate(phenotype_pos = as.integer(ceiling((gene_end + gene_start)/2))) %>%
    dplyr::select(required_gene_meta_columns, everything()) %>% #Reorder columns
    as.data.frame()
  rownames(gene_data) = gene_data$gene_id

  #identify shared genes
  shared_genes = intersect(featureCounts$gene_id, gene_data$gene_id)
  gene_data = gene_data[shared_genes,]

  #Make a read counts matrix
  count_matrix = dplyr::select(featureCounts, -gene_id, length)
  count_matrix = as.matrix(count_matrix)
  row.names(count_matrix) = featureCounts$gene_id
  count_matrix = count_matrix[shared_genes,]

  #Add gene lengths to the gene metadata file
  length_df = dplyr::transmute(read_counts, gene_id, phenotype_length = length)
  gene_data = dplyr::left_join(gene_data, length_df, by = "gene_id") %>%
    as.data.frame()
  rownames(gene_data) = gene_data$gene_id

  #Identify shared samples
  shared_samples = intersect(sample_metadata$sample_id, colnames(count_matrix))
  count_matrix = count_matrix[,shared_samples]

  #Prepare sample metadata
  sample_meta = as.data.frame(sample_metadata)
  row.names(sample_meta) = sample_meta$sample_id
  sample_meta = sample_meta[shared_samples,]

  #Make a summarizedExperiment object
  se = SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = count_matrix),
    colData = sample_meta,
    rowData = gene_data)

  return(se)
}

makeSummarizedExperiment <- function(assay, row_data, col_data, assay_name){

  #Make dfs
  row_df = as.data.frame(row_data)
  rownames(row_df) = row_data$phenotype_id

  col_df = as.data.frame(col_data)
  rownames(col_df) = col_data$sample_id

  #Keep only overlapping samples
  shared_samples <- BiocGenerics::intersect(col_df$sample_id, colnames(assay))
  col_df <- col_df[shared_samples,]
  assay = assay[,shared_samples]

  #Make assay list
  assay = as.matrix(assay)
  assay_list = list()
  #Reorder columns and rows to match metadata
  assay_list[[assay_name]] = assay[rownames(row_df), rownames(col_df)]


  #Make a summarizedExperiment object
  se = SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    colData = col_df,
    rowData = row_df)
}

#' Make Summarized Experiment object from raw count matrix
#'
#' @param assay count matrix of any phenotype quantification.
#' @param row_data phenotype metadata
#' @param col_data sample metadata
#' @param assay_name custom assay name to be created in SE. Default "counts"
#' @param quant_method Quantification method. Can be gene_counts, leafcutter, txrevise, transcript_usage or exon_counts
#' @author Nurlan Kerimov
#' @export
makeSummarizedExperimentFromCountMatrix <- function(assay, row_data, col_data, assay_name = "counts", quant_method = "gene_counts"){
  #Make dfs
  row_df = BiocGenerics::as.data.frame(row_data)
  rownames(row_df) = row_data$phenotype_id

  col_df = BiocGenerics::as.data.frame(col_data)
  rownames(col_df) = col_data$sample_id

  shared_samples <- BiocGenerics::intersect(col_df$sample_id, colnames(assay))
  col_df <- col_df[shared_samples,]

  assay = assay %>% dplyr::filter(!(phenotype_id %like% "PAR_Y")) %>% reformatPhenotypeId(quant_method = quant_method)
  rownames(assay) <- assay$phenotype_id
  assay = assay[,shared_samples]

  if (quant_method %in% c("exon_counts","transcript_usage", "txrevise")) {
    # remove invalid gene types from assay
    valid_gene_types = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                         "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                         "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                         "antisense","sense_intronic","sense_overlapping")
    row_df = row_df[row_df$gene_type %in% valid_gene_types,]

    if (quant_method == "exon_counts") {
      # remove the exons which only less than 5 samples have overlapped read(s)
      assay = assay[rowSums(assay) > 5,]
    }

    shared_phenotypes <- BiocGenerics::intersect(rownames(row_df), rownames(assay))
    assay = assay[shared_phenotypes,]
    row_df = row_df[shared_phenotypes,]

    if (quant_method %in% c("transcript_usage","txrevise")) {
      assay_name <- "tpms"
    }
  }

  dummy <- assertthat::assert_that(all(rownames(assay) %in% row_df$phenotype_id), msg = "Some phenotypes in assay missing metadata information")

  #Make assay list
  assay = as.matrix(assay)
  assay_list = list()
  #Reorder columns and rows to match metadata
  assay_list[[assay_name]] = assay[rownames(row_df), rownames(col_df)]

  #Make a summarizedExperiment object
  se = SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    colData = col_df,
    rowData = row_df)
}

#' Filter a Summarized Experiment according to provided parameters
#'
#' @param se Summarized Experiment object to be filteres
#' @param valid_chromosomes list (string vector) of valid gene chromosomes
#' @param valid_gene_types list (string vector) of valid gene types
#' @param filter_rna_qc Boolean value
#' @param filter_genotype_qc Boolean value
#' @author Kaur Alasoo
#' @export
filterSummarizedExperiment <- function(se, valid_chromosomes = NA, valid_gene_types = NA, filter_rna_qc = FALSE, filter_genotype_qc = FALSE){

  #Filter chromosomes
  if(!is.na(valid_chromosomes[1])){
    assertthat::assert_that(assertthat::has_name(SummarizedExperiment::rowData(se), "chromosome"))
    se = se[SummarizedExperiment::rowData(se)$chromosome %in% valid_chromosomes,]
  }
  #Filter gene types
  if(!is.na(valid_gene_types[1])){
    assertthat::assert_that(assertthat::has_name(SummarizedExperiment::rowData(se), "gene_type"))
    se = se[SummarizedExperiment::rowData(se)$gene_type %in% valid_gene_types,]
  }

  #Filter RNA QC
  if(filter_rna_qc){
    se = se[,se$rna_qc_passed]
  }

  #Filter genotype QC
  if(filter_genotype_qc){
    se = se[,se$genotype_qc_passed]
  }
  return(se)
}

#' Normalises the assay of SE with cqn method
#'
#' @param se Summarized Experiment object to
#' @param assay_name custom assay_name (default: counts)
#' @author Kaur Alasoo
#' @export
normaliseSE_cqn <- function(se, assay_name = "counts"){
  require("cqn")

  #Extract fields
  row_data = SummarizedExperiment::rowData(se)
  col_data = SummarizedExperiment::colData(se)
  assays = SummarizedExperiment::assays(se)

  #Ensure that required columns are present
  assertthat::assert_that(assertthat::has_name(row_data, "phenotype_length"))
  assertthat::assert_that(assertthat::has_name(row_data, "phenotype_gc_content"))

  #Perform cqn-normalisation
  count_matrix = assays[[assay_name]]
  expression_cqn = cqn::cqn(counts = count_matrix,
                            x = row_data$phenotype_gc_content,
                            lengths = row_data$phenotype_length, verbose = TRUE)
  cqn_matrix = expression_cqn$y + expression_cqn$offset

  #Update assays
  assays[["cqn"]] = cqn_matrix

  #Make ab update se object
  se = SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = col_data,
    rowData = row_data)
}

#' Normalises the assay of SE with TPM method
#'
#' @param se Summarized Experiment object
#' @param assay_name custom assay_name (default: counts)
#' @param fragment_length custom fragment_length (default: 250)
#' @author Kaur Alasoo
#' @export
normaliseSE_tpm <- function(se, assay_name = "counts", fragment_length = 250){

  #Extract fields
  row_data = SummarizedExperiment::rowData(se)
  col_data = SummarizedExperiment::colData(se)
  assays = SummarizedExperiment::assays(se)

  #Ensure that required columns are present
  assertthat::assert_that(assertthat::has_name(row_data, "phenotype_length"))

  #Perform tpm-normalisation
  count_matrix = assays[[assay_name]]
  lengths = row_data$phenotype_length
  scaling_factor = colSums((count_matrix*fragment_length)/lengths)

  tpm = t(t(count_matrix * fragment_length * 1e6)/scaling_factor)
  tpm = tpm/lengths

  #Update assays
  assays[["tpms"]] = tpm

  #Make ab update se object
  se = SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = col_data,
    rowData = row_data)
}

filterSE_expressedGenes <- function(se, min_count = 1, min_fraction = 0.1, assay_name = "counts"){

  #Count the number of samples with at least min_count reads for each gene
  count_matrix = assays(se)[[assay_name]]
  expressed_counts = rowSums(count_matrix >= min_count)
  thresh_count = ceiling(0.1*(dim(count_matrix)[2]))
  expressed_genes = names(which(expressed_counts > thresh_count))

  #Filter the se
  final_se = se[expressed_genes,]
  return(final_se)
}

#' Quantile normalise SummarizedExperiment by rows
#'
#' @param se SummarizedExperiment object
#' @param assay_name Name of the assay to be processed in the se object
#' @param n_pcs Number of desired principal components
#' @param log_transform Boolean value
#' @param column_prefix custom column_prefix
#' @param feature_id default: sample_id
#'
#' @return SummarizedExperiment object with quantile-normalised data in the qnorm assay
#' @export
transformSE_PCA <- function(se, assay_name = "cqn", n_pcs = NULL, log_transform = FALSE, column_prefix = "", feature_id = "sample_id", ...){
  #Extract data
  matrix = SummarizedExperiment::assays(se)[[assay_name]]
  sample_data = SummarizedExperiment::colData(se) %>% as.data.frame() %>% dplyr::as_tibble()

  if(log_transform){
    matrix = log(matrix + 0.1, 2)
  }

  #Perform PCA of gene expression matrix add experimental design metadata to the results
  pca = prcomp(t(matrix), ...)
  if(is.null(n_pcs)){
    n_pcs = ncol(matrix)
  }
  pca_mat = as.data.frame(pca$x[,1:n_pcs])
  colnames(pca_mat) = paste0(column_prefix, colnames(pca_mat))
  pca_matrix = pca_mat %>%
    dplyr::mutate(sample_id = rownames(pca$x)) %>%
    dplyr::rename_(.dots = setNames("sample_id", feature_id)) %>% #Hack to make renaming work
    dplyr::left_join(sample_data, by = feature_id)
  #Calculate variance explained by each component
  var_exp = (pca$sdev^2) / sum(pca$sdev^2)
  return(list(pca_matrix = pca_matrix, pca_object = pca, var_exp = var_exp))
}


removeGeneVersion <- function(read_counts){
  read_counts$gene_id = (dplyr::select(read_counts, gene_id) %>% tidyr::separate(gene_id, c("gene_id", "suffix"), sep = "\\."))$gene_id
  return(read_counts)
}

#' Reformat phenotype IDs according to quantification method
#'
#' @param read_counts count matrix of any phenotype quantification.
#' @param quant_method Quantification method. Can be gene_counts, leafcutter, txrevise, transcript_usage or exon_counts
#' @author Nurlan Kerimov
#' @export
reformatPhenotypeId <- function(read_counts, quant_method = "gene_counts") {
  if (quant_method %in% c("gene_counts","transcript_usage")) {
    read_counts$phenotype_id = (dplyr::select(read_counts, phenotype_id) %>% tidyr::separate(phenotype_id, c("phenotype_id", "suffix"), sep = "\\."))$phenotype_id
  }
  return(read_counts)
}

mergeCountsSEs <- function(se1, se2){
  #Extract fields
  counts1 = SummarizedExperiment::assays(se1)$counts
  col1 = SummarizedExperiment::colData(se1)

  counts2 = SummarizedExperiment::assays(se2)$counts
  col2 = SummarizedExperiment::colData(se2)

  #Keep only shared gene ids
  shared_genes = intersect(rownames(counts1), rownames(counts2))

  #Merge data
  shared_cols = intersect(colnames(col1), colnames(col2))
  merged_coldata = rbind(col1[,shared_cols], col2[,shared_cols])
  merged_counts = cbind(counts1[shared_genes,], counts2[shared_genes,])

  #Extraxct row dara from first SE
  row_data = SummarizedExperiment::rowData(se1[shared_genes,])

  se = SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = merged_counts),
    colData = merged_coldata,
    rowData = row_data)
  return(se)
}

#' Extract a subset of the data from a SummarizedExperiment based on column information
#'
#' @param se SummarizedExperiment object
#' @param column Name of the column in the colData object
#' @param value Value of the column to be extracted
#'
#' @return Subsetted SummmarizedExperiment object
#' @export
subsetSEByColumnValue <- function(se, column, value){
  selection = SummarizedExperiment::colData(se)[,column] == value
  result = se[,selection]
  return(result)
}

#' normaliseSE_ratios
#'
#' @param se SummarizedExperiment object
#' @param assay_name custom assay name for generated assay
#'
#' @return SummarizedExperiment object with additional assay
#' @export
normaliseSE_ratios <- function(se, assay_name = "tpms"){

  #Extract rowData and check for required columns
  row_data = SummarizedExperiment::rowData(se)
  assertthat::assert_that(assertthat::has_name(row_data, "phenotype_id"))
  assertthat::assert_that(assertthat::has_name(row_data, "quant_id"))

  #Extract metadata
  tx_map = row_data %>%
    as.data.frame() %>%
    dplyr::as_tibble() %>%
    dplyr::select(phenotype_id, quant_id)

  #Extract assays
  assay_list = SummarizedExperiment::assays(se)
  assay_matrix = assay_list[[assay_name]]

  #Calculate transcript ratios
  transcript_ratios = calculateTranscriptUsage(assay_matrix, tx_map) %>%
    replaceNAsWithRowMeans()

  #Reorder rows
  reordered_ratios = transcript_ratios[row_data$phenotype_id,]

  #Update assays
  assay_list[["usage"]] = reordered_ratios

  #Make ab update se object
  se = SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    colData = SummarizedExperiment::colData(se),
    rowData = row_data)
  return(se)
}


#' Quantile normalise SummarizedExperiment by rows
#'
#' @param se SummarizedExperiment object
#' @param assay_name Name of the assay to be normalised in the se object
#'
#' @return SummarizedExperiment object with quantile-normalised data in the qnorm assay
#' @export
normaliseSE_quantile <- function(se, assay_name = "usage"){

  #Extract assays
  assay_list = SummarizedExperiment::assays(se)
  assay_matrix = assay_list[[assay_name]]

  #Quantile normalise
  qnorm = quantileNormaliseRows(assay_matrix)
  assay_list[["qnorm"]] = qnorm

  #Make ab update se object
  se = SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    colData = SummarizedExperiment::colData(se),
    rowData = SummarizedExperiment::rowData(se))
  return(se)
}

#' Remove phenotypes from SummarizedExperiments that do not have cis genotypes
#'
#' @param se SummarizedExperiment object
#' @param variant_information Variant information data frame from importVariantInformation()
#' @param variant_information Maximal distance from phenotype_pos
#' @param variant_information Minimal number of variants in cis before the phenotype is removed.
#'
#' @return filtered SummarizedExperiment object where phenotypes with no cis variants are removed.
#' @export
checkCisVariants <- function(se,
                             variant_information,
                             cis_distance = 1000000,
                             min_cis_variant = 5) {

  #Extract gene metadata from SummarizedExperiment
  gene_data = SummarizedExperiment::rowData(se) %>% as.data.frame() %>% as_tibble()

  #Check that all required columns are there
  assertthat::assert_that(assertthat::has_name(gene_data, "chromosome"))
  assertthat::assert_that(assertthat::has_name(gene_data, "phenotype_pos"))
  assertthat::assert_that(assertthat::has_name(gene_data, "strand"))

  assertthat::assert_that(assertthat::has_name(variant_information, "chr"))
  assertthat::assert_that(assertthat::has_name(variant_information, "pos"))

  #Make GRanges objects
  gene_ranges = GenomicRanges::GRanges(
    seqnames = gene_data$chromosome,
    ranges = IRanges::IRanges(
      start = gene_data$phenotype_pos - cis_distance,
      end = gene_data$phenotype_pos + cis_distance
    ),
    strand = "*"
  )
  var_ranges = GenomicRanges::GRanges(
    seqnames = var_info$chr,
    ranges = IRanges::IRanges(start = var_info$pos, end = var_info$pos),
    strand = "+"
  )
  olap_count = GenomicRanges::countOverlaps(gene_ranges, var_ranges, ignore.strand = TRUE)
  count_df = dplyr::select(gene_data, phenotype_id, chromosome, phenotype_pos) %>%
    dplyr::mutate(snp_count = olap_count) %>%
    dplyr::arrange(snp_count) %>%
    dplyr::filter(snp_count >= min_cis_variant)

  #Keep phenotypes that contain enough variants nearby
  se = se[count_df$phenotype_id,]
  return(se)
}

extractPhentypeFromSE <- function(phenotype_id, se, assay){

  #extract single phenotype
  pheno_mat = assays(se)[[assay]]
  pheno_row = pheno_mat[phenotype_id,]
  pheno_df = dplyr::data_frame(phenotype_id = phenotype_id, sample_id = names(pheno_row), phenotype_value = pheno_row)

  #Extract sample metadata
  sample_meta = SummarizedExperiment::colData(se) %>% as.data.frame() %>% dplyr::as_tibble()
  pheno_df = dplyr::left_join(pheno_df, sample_meta, by = "sample_id")
  return(pheno_df)
}
