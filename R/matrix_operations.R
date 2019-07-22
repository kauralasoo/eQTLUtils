calculateMean <- function(matrix, design, factor, sample_id_col = "sample_id", na.rm = FALSE){
  #Calculate the mean value in matrix over all possible factor values.

  #If the factor is not a factor then make it a factor.
  if(!is.factor(design[,factor])){
    design[,factor] = factor(design[,factor])
  }

  #Set sample_id column as rownames
  rownames(design) = design[,sample_id_col]
  factor = design[,factor]
  levs = levels(factor)
  result = c()
  for (lev in levs){
    filter = factor == lev
    samples = rownames(design[filter,])
    mat = matrix[,samples]
    mat = rowMeans(mat, na.rm)
    result = cbind(result, mat)
  }
  colnames(result) = levs
  return(data.frame(result))
}

zScoreNormalize <- function(matrix){
  #Normalize expression matrix by z-score
  matrix = matrix - rowMeans(matrix)
  matrix = matrix / apply(matrix, 1, sd)
  return(matrix)
}

#' replaceNAsWithRowMeans
#'
#' @param matrix matrix
#'
#' @export
replaceNAsWithRowMeans <- function(matrix){
  #replace with row means
  na_pos = which(is.na(matrix), arr.ind = TRUE)
  matrix[na_pos] = rowMeans(matrix, na.rm=TRUE)[na_pos[,1]]

  #If there are addional NAs left (whole row NAs) then replace with 0
  matrix[is.na(matrix)] = 0
  return(matrix)
}

#' Force a vector of values into standard normal distribution
#'
#' @param x numeric vector with arbitrary distribution
#'
#' @return Vector with a standard normal distribution
#' @export
quantileNormaliseVector = function(x){
  qnorm(rank(x,ties.method = "random")/(length(x)+1))
}

quantileNormaliseMatrix <- function(matrix){
  quantile_matrix = matrix(0, nrow(matrix), ncol(matrix))
  for (i in seq_along(matrix[1,])){
    quantile_matrix[,i] = quantileNormaliseVector(matrix[,i])
  }
  #Add names
  rownames(quantile_matrix) = rownames(matrix)
  colnames(quantile_matrix) = colnames(matrix)
  return(quantile_matrix)
}

quantileNormaliseCols <- function(matrix,...){
  quantileNormaliseMatrix(matrix, ...)
}

quantileNormaliseRows <- function(matrix,...){
  t(quantileNormaliseMatrix(t(matrix), ...))
}

#' Calculate the transcript usage
#'
#' @param expression_matrix expression_matrix 
#' @param phenotype_map phenotype_map 
#'
#' @export
calculateTranscriptUsage <- function(expression_matrix, phenotype_map){

  #Check that gene_transcript_map has the correct columns
  assertthat::assert_that(assertthat::has_name(phenotype_map, "quant_id"))
  assertthat::assert_that(assertthat::has_name(phenotype_map, "phenotype_id"))
  assertthat::assert_that(length(colnames(phenotype_map)) == 2)

  #group phenotype_ids by quant_id
  groups = dplyr::group_by(phenotype_map, quant_id) %>% tidyr::nest()

  #Extract individual phenotypes
  matrix_list = purrr::map(groups$data, ~expression_matrix[.$phenotype_id,,drop = FALSE])

  #Divide by row sums
  usage_list = purrr::map(matrix_list, ~t(t(.)/apply(.,2, sum)))
  transcript_usage = do.call(rbind, usage_list)

  return(transcript_usage)
}
