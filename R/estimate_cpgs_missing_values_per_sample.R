#' Estimate the proportion of CpG sites that are missing values in each sample
#' 
#' @param meth_rse A RangedSummarizedExperiment with methylation values.
#' @param assay The assay from meth_rse to extract values from. 
#' Should be either an index or the name of an assay. Default is the first assay. 
#' @param number_cpgs Number of CpG sites to sample. Default is 1000. 
#' @param first_seq_only A logical value indicating whether to restrict sampling to 
#' the first sequence level of meth_rse to increase speed. Default value is TRUE. 
#' @return A numeric vector with the proportion of missing values found in each sample.
#' @export
estimate_cpgs_missing_values_per_sample = function(meth_rse, assay = 1, number_cpgs = 1000, first_seq_only = TRUE){
  
  # Extract CpGs from meth_rse as a GRanges
  cpgs = rowRanges(meth_rse)
  
  # If first_seq_only is set to TRUE, subset cpgs for ranges on the first sequence
  if(first_seq_only){
    cpgs = cpgs[seqnames(cpgs) == seqlevels(cpgs)[1]]
  }
  
  # Sample CpGs
  sample_cpgs = sample(cpgs, number_cpgs)
  
  # Subset meth_rse for the sampled CpGs
  meth_rse_sample = subsetByOverlaps(meth_rse, sample_cpgs)
  
  # Calculate the proportion of missing values per sample and return
  missing_values = colSums(is.na(assay(meth_rse_sample, i = assay)))/number_cpgs
  return(missing_values)
  
}
