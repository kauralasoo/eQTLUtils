% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qtl_QTLTools.R
\name{saveQTLToolsMatrices}
\alias{saveQTLToolsMatrices}
\title{Save a list of matrices into a suitable format for QTLTools}
\usage{
saveQTLToolsMatrices(
  data_list,
  output_dir,
  file_suffix = "bed",
  file_prefix = "",
  col_names = TRUE
)
}
\arguments{
\item{data_list}{list of matrices}

\item{output_dir}{relative path to the output dir}

\item{file_suffix}{suffix added to each file after their name in the list.}
}
\value{
None
}
\description{
Works with expression and covariates matrices.
}
\author{
Kaur Alasoo
}
