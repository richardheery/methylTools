% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_cpgs_missing_values_per_sample.R
\name{estimate_cpgs_missing_values_per_sample}
\alias{estimate_cpgs_missing_values_per_sample}
\title{Estimate the proportion of CpG sites that are missing values in each sample}
\usage{
estimate_cpgs_missing_values_per_sample(
  meth_rse,
  assay = 1,
  number_cpgs = 1000,
  filter_seq_levels = NULL
)
}
\arguments{
\item{meth_rse}{A RangedSummarizedExperiment with methylation values.}

\item{assay}{The assay from meth_rse to extract values from.
Should be either an index or the name of an assay. Default is the first assay.}

\item{number_cpgs}{Number of CpG sites to sample. Default is 1000.}

\item{filter_seq_levels}{An optional vector with the sequence levels to restrict sampling of CpGs to.
Setting to just a single sequence level can substantially speed up running time.}
}
\value{
A numeric vector with the proportion of missing values found in each sample.
}
\description{
Estimate the proportion of CpG sites that are missing values in each sample
}
