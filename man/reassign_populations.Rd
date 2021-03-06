% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pop_assign.R
\name{reassign_populations}
\alias{reassign_populations}
\title{Reassigns already assigned populations according to new thresholds}
\usage{
reassign_populations(distance_matrix, abs_threshold, rel_threshold)
}
\arguments{
\item{distance_matrix}{distance matrix dataframe created by assign_populations() function}

\item{abs_threshold}{sample's minimum distance to any reference population should be less than this threshold to be assigned to population}

\item{rel_threshold}{sample's minimum distance to any reference population should be rel_threshold times smaller than second minimum distanceto be assigned to population}
}
\description{
Reassigns already assigned populations according to new thresholds
}
\author{
Nurlan Kerimov
}
