qtltools_prepareFeatureCountsData <- function(se){
  require("cqn")
  
  #Specify valid chromsomes and valid gene types
  valid_gene_types = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                       "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                       "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                       "antisense","sense_intronic","sense_overlapping")
  valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                        "2","20","21","22","3","4","5","6","7","8","9")
  
  #Filter SE to keep correct chromosomes and QC-passed samples
  se_filtered = filterSummarizedExperiment(se, valid_chromosomes = valid_chromosomes,
                                           valid_gene_types = valid_gene_types,
                                           filter_rna_qc = TRUE, filter_genotype_qc = TRUE)
  
  #Identify expressed genes (At least 1 count in 10% of the samples)
  se_expressed = filterSE_expressedGenes(se_filtered, min_count = 1, min_fraction = 0.1)
  
  #Normalise and make QTLtools matrix
  cqn_se = normaliseSE_cqn(se_expressed)
  return(cqn_se)
}