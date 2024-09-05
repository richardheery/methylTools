#' Estimate the proportion of CpG sites that are missing values in each sample
#' 
#' @param meth_rse A RangedSummarizedExperiment with methylation values.
#' @param assay The assay from meth_rse to extract values from. 
#' Should be either an index or the name of an assay. Default is the first assay. 
#' @param number_cpgs Number of CpG sites to sample. Default is 1000. 
#' @param filter_seq_levels An optional vector with the sequence levels to restrict sampling of CpGs to.
#' Setting to just a single sequence level can substantially speed up running time.
#' @return A numeric vector with the proportion of missing values found in each sample.
#' @export
estimate_cpgs_missing_values_per_sample = function(meth_rse, assay = 1, number_cpgs = 1000, filter_seq_levels = NULL){
  
  # Extract CpGs from meth_rse as a GRanges
  cpgs = rowRanges(meth_rse)
  
  # If filter_seq_levels is NULL, set to all seqlevels from meth_rse
  if(is.null(filter_seq_levels)){
    filter_seq_levels = seqlevels(meth_rse)
  }
  
  # Subset CpGs for those on filter_seq_levels
  cpgs = cpgs[seqnames(cpgs) %in% filter_seq_levels]
  
  # Sample CpGs
  sample_cpgs = sample(cpgs, number_cpgs)
  
  # Subset meth_rse for the sampled CpGs
  meth_rse_sample = subsetByOverlaps(meth_rse, sample_cpgs)
  
  # Calculate the proportion of missing values per sample and return
  missing_values = MatrixGenerics::colSums2(is.na(assay(meth_rse_sample, i = assay)))/number_cpgs
  return(missing_values)
  
}
