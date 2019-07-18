#' Import variant information from GDS file
#'
#' @param gdsfile Path to the genotype file in GDS format.
#'
#' @return data frame with variant information
#' @export
importVariantInformationFromGDS <- function(gdsfile){

  #Import individual columns
  snp_pos = GDSArray::GDSArray(gdsfile, "snp.position")
  snp_chromosome = GDSArray::GDSArray(gdsfile, "snp.chromosome")
  snp_id = GDSArray::GDSArray(gdsfile, "snp.rs.id")

  #Make a data frame
  snp_df = dplyr::data_frame(gds_snp_id = as.integer(names(snp_id)),
                             chromosome = as.vector(snp_chromosome),
                             pos = as.vector(snp_pos),
                             snp_id = as.vector(snp_id))
  return(snp_df)
}

#' Extract genotypes of a single genetic variant from the GDS file
#'
#' @param variant_id Id of the genetic variant
#' @param variant_information Variant infromation constructed by importVariantInformationFromGDS
#' @param gdsfile Path to the GDS file
#'
#' @return data frame with the genotypes
#' @export
extractVariantGenotypeFromGDS <- function(variant_id, variant_information, gdsfile){
  assertthat::assert_that(length(variant_id) == 1)

  #Extract variant gds_id from variant infromation
  var_filtered = dplyr::filter(variant_information, snp_id == variant_id)
  selected_gds_id = var_filtered$gds_snp_id
  assertthat::assert_that(length(selected_gds_id) == 1)

  #Extract genotype from the gds file
  geno = GDSArray::GDSArray(gdsfile, "genotype")
  genotype = geno[selected_gds_id,]
  genotype_df = dplyr::data_frame(genotype_id = colnames(geno), genotype_value = genotype, snp_id = variant_id)
  return(genotype_df)
}

#' Extract a genotype matrix corresponding to a specific genomic region
#'
#' @param chr Chromosome of the region
#' @param start Start coordinate of the region
#' @param end End coordinate of the reion
#' @param variant_information Variant infromation constructed by importVariantInformationFromGDS
#' @param gdsfile Path to the GDS file
#'
#' @return genotype matrix
#' @export
extractGenotypeMatrixFromGDS <- function(chr, start, end, variant_information, gdsfile){

  #Extract variant ids from variant infromation
  var_meta = dplyr::filter(variant_information, chromosome == chr, pos > start, pos < end)
  gds_ids = var_meta$gds_snp_id

  #Extract genotype from the gds file
  geno = GDSArray::GDSArray(gdsfile, "genotype")
  genotype = as.matrix(geno[gds_ids,])
  rownames(genotype) = var_meta$snp_id

  return(genotype)
}
