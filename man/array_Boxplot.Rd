% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/array_qc_utils.R
\name{array_Boxplot}
\alias{array_Boxplot}
\title{Generate Boxplot QC Plot for n random samples}
\usage{
array_Boxplot(
  se,
  n = 10,
  assay_name = "exprs",
  filter_quality = TRUE,
  pdf_output = FALSE,
  output_dir = "./"
)
}
\arguments{
\item{se}{SummarizedExperiment file to be analysed}

\item{n}{Number of random samples (Default:10)}

\item{assay_name}{Assay name to be analysed (exprs, norm_exprs or batch_exprs)}

\item{filter_quality}{Boolean value if samples with RNA_QC_passed = FALSE should be excluded (Default:TRUE)}

\item{pdf_output}{Boolean value if pdf output should be created (Default:FALSE)}

\item{output_dir}{html file output dir, if html_output is TRUE (Default:current directory)}
}
\value{
Boxplot plot of sample expression values
}
\description{
Generate Boxplot QC Plot for n random samples
}
\author{
Liis Kolberg
}
