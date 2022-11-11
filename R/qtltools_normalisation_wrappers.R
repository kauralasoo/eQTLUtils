#' Normalise phenotype SummarisedExperiment object for QTL mapping with QTLtools
#'
#' @param se SummarizedExperiment object used for QTL mapping
#' @param quant_method Quantification method used to generate the data. Valid options are: featureCounts, array, leafcutter
#'
#' @return Normalised SummarizedExperiment object
#' @export
qtltoolsPrepareSE <- function(se, quant_method, filter_rna_qc = TRUE, filter_genotype_qc = TRUE, keep_XY = FALSE){
  require("cqn")

  #Check for valid quant_methods
  valid_quant_methods = c("gene_counts", "exon_counts","HumanHT-12_V4", "leafcutter", "txrevise", "transcript_usage")
  assertthat::assert_that(quant_method %in% valid_quant_methods, msg = "The specified quant_method is not valid.")

  #Specify valid chromsomes and valid gene types
  valid_gene_types = c("lincRNA","lncRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                       "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                       "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                       "antisense","sense_intronic","sense_overlapping")
  valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                        "2","20","21","22","3","4","5","6","7","8","9")
  if(keep_XY){
    valid_chromosomes = c(valid_chromosomes, c("X","Y"))
  }

  #Use different normalisation strategy for each quantification method
  if(quant_method %in% c("exon_counts", "gene_counts")){
    #Normalise featureCounts data

    #Filter SE to keep correct chromosomes and QC-passed samples
    se_filtered = eQTLUtils::filterSummarizedExperiment(se, valid_chromosomes = valid_chromosomes,
                                             valid_gene_types = valid_gene_types,
                                             filter_rna_qc = filter_rna_qc, filter_genotype_qc = filter_genotype_qc)

    #Normalise and make QTLtools matrix
    se_norm = normaliseSE_cqn(se_filtered, assay_name = "counts")

    } else if(quant_method == "HumanHT-12_V4"){

      #Filter SE to keep correct chromosomes and QC-passed samples
      message(" ## Filtering out invalid chromosomes, RNA_QC failed samples and genotype_QC failed samples ")
      se_filtered = eQTLUtils::filterSummarizedExperiment(se, valid_chromosomes = valid_chromosomes,
                                                          valid_gene_types = valid_gene_types,
                                                          filter_rna_qc = filter_rna_qc, filter_genotype_qc = filter_genotype_qc)

      #Normalize and regress out batch effects
      message(" ## Normalizing and regressing out the batch effects")
      se_norm = eQTLUtils::array_normaliseSE2(se_filtered, norm_method = "quantile", assay_name = "exprs",
                                           log_transform = TRUE, adjust_batch = TRUE)

    } else if(quant_method == "leafcutter"){
      #Normalise LeafCutter exon exicision count

      #Keep events that did not map to any gene
      valid_gene_types = c(valid_gene_types, "leafcutter")

      #Filter se to keep correct chromosomes and gene biotypes
      se_filtered = eQTLUtils::filterSummarizedExperiment(se, valid_chromosomes = valid_chromosomes,
                                                          valid_gene_types = valid_gene_types,
                                                          filter_rna_qc = filter_rna_qc, filter_genotype_qc = filter_genotype_qc)

      message("Calculating ratios...")
      se_ratios = eQTLUtils::normaliseSE_ratios(se_filtered, assay_name = "counts")

      message("Performing inverse normal transformation...")
      se_norm = eQTLUtils::normaliseSE_quantile(se_ratios, assay_name = "usage")

    } else if(quant_method %in% c("txrevise", "transcript_usage")){

      #Filter SE to keep correct chromosomes and QC-passed samples
      se_filtered = eQTLUtils::filterSummarizedExperiment(se, valid_chromosomes = valid_chromosomes,
                                               valid_gene_types = valid_gene_types,
                                               filter_rna_qc = filter_rna_qc, filter_genotype_qc = filter_genotype_qc)

      message("Calculating ratios...")
      se_ratios = eQTLUtils::normaliseSE_ratios(se_filtered, assay_name = "tpms")

      message("Performing inverse normal transformation...")
      se_norm = eQTLUtils::normaliseSE_quantile(se_ratios, assay_name = "usage")

    }

  return(se_norm)
}
