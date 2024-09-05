#' Add number of methylated and unmethylated read counts to a methylation RSE with methylation proportion and coverage
#' 
#' @param meth_rse A RangedSummarizedExperiment with methylation values.
#' @param proportion_assay The assay of meth_rse which corresponds to the proportion of methylated reads. 
#' Can be either a numeric index or the name of the assay. Default is the first assay.
#' @param coverage_assay The assay of meth_rse which corresponds to the coverage (the total number of reads).
#' Can be either a numeric index or the name of the assay. Default is the second assay.
#' @return A RangedSummarized experiment identical to meth_rse with two additional assays added: one for methylated reads and another for unmethylated reads.
#' @export
add_counts_to_meth_rse = function(meth_rse, proportion_assay = 1, coverage_assay = 2){
  
  # Add assay with the number of methylated reads (the coverage multiplied by the proportion of methylated reads)
  assays(meth_rse)[["meth_reads"]] = assays(meth_rse)[[coverage_assay]] * assays(meth_rse)[[proportion_assay]]
  
  # Add asay with number of unmethylated reads (the coverage less the number of methylated reads)
  assays(meth_rse)[["unmeth_reads"]] =  assays(meth_rse)[[coverage_assay]] - assays(meth_rse)[["meth_reads"]]
  
  # Return meth_rse
  return(meth_rse)
  
}

#' Test differential methylation in predefined genomic regions using methylSig
#' 
#' First calculates the total number of methylated reads and the total number of reads
#' for all CpG sites in each region in each sample before testing differential methylation using 
#' the beta-binomial test approach from methylSig.  
#' 
#' @param meth_rse A RangedSummarizedExperiment with methylation values.
#' @param genomic_regions GRanges object with regions to get methylated and unmethylated counts for.  
#' @param meth_reads_assay The assay of meth_rse which corresponds to the read counts for methylated reads. 
#' Can be either a numeric index or the name of the assay. Default is "meth_reads".
#' @param coverage_assay The assay of meth_rse which corresponds to the total read counts for methylated and unmethylated reads. 
#' @param max_sites_per_chunk The approximate maximum number of methylation sites to try to load into memory at once. 
#' The actual number loaded may vary depending on the number of methylation sites overlapping each region, 
#' but so long as the size of any individual regions is not enormous (>= several MB), it should vary only very slightly. 
#' Some experimentation may be needed to choose an optimal value as low values will result in increased running time, 
#' while high values will result in a large memory footprint without much improvement in running time. 
#' Default is floor(62500000/ncol(meth_rse)), resulting in each chunk requiring approximately 500 MB of RAM. 
#' @param na.rm TRUE or FALSE indicating whether to remove NA values when calculating summaries. Default value is TRUE. 
#' @param group_column Name of column in colData(meth_rse) indicating the groups to be compared. 
#' @param case Which level of group_column corresponds to the group of interest.
#' @param control Which level of group_column corresponds to the control group. 
#' @param BPPARAM A BiocParallelParam object. Defaults to `BiocParallel::bpparam()`. 
#' @return A GRanges object with the results of methylSig::diff_methylsig(). 
#' @export
diff_meth_methylsig = function(meth_rse, genomic_regions, meth_reads_assay = "meth_reads", coverage_assay = "cov", 
  max_sites_per_chunk = NULL, na.rm = TRUE, group_column, case, control, BPPARAM = BiocParallel::bpparam()){
  
  # Check that group_column is in colData(meth_Rse) and that it contains the levels indicated by case and control
  if(!group_column %in% names(SummarizedExperiment::colData(meth_rse))){
    stop(paste(group_column, "is not a column in colData(meth_rse)"))
  } else {
    if(!case %in% SummarizedExperiment::colData(meth_rse)[[group_column]] || !control %in% SummarizedExperiment::colData(meth_rse)[[group_column]]){
      stop(paste("Either", case, "or", control, "is not present in meth_rse"))
    }
  }
  
  # Check that meth_reads_assay and coverage_assay present in meth_rse
  if(!any(c(meth_reads_assay, coverage_assay) %in% names(SummarizedExperiment::assays(meth_rse)))){
    stop(paste("Either", meth_reads_assay, "or", coverage_assay, "not present in names(assays(meth_rse))"))
  }
  
  # Get methylated reads for genomic regions
  message("Summing methylated reads for genomic regions")
  genomic_region_meth_reads = methodical::summarizeRegionMethylation(meth_rse = meth_rse, assay = meth_reads_assay, 
    genomic_regions = genomic_regions, col_summary_function = "colSums2", 
    max_sites_per_chunk = max_sites_per_chunk, na.rm = na.rm, BPPARAM = BPPARAM)
  
  # Get Coverage for genomic regions
  message("Summing coverage for genomic regions")
  genomic_region_coverage = methodical::summarizeRegionMethylation(meth_rse = meth_rse, assay = coverage_assay, 
    genomic_regions = genomic_regions, col_summary_function = "colSums2", 
    max_sites_per_chunk = max_sites_per_chunk, na.rm = na.rm, BPPARAM = BPPARAM)
  
  # Set NA values equal to 0 
  genomic_region_meth_reads[is.na(genomic_region_meth_reads)] = 0
  genomic_region_coverage[is.na(genomic_region_coverage)] = 0
  
  # Create a BSSeq object
  bsseq = bsseq::BSseq(M = as.matrix(genomic_region_meth_reads), 
    Cov = as.matrix(genomic_region_coverage), gr = genomic_regions)
  pData(bsseq) = colData(meth_rse)
  
  # Test differential methylation and return
  diff_meth_results = methylSig::diff_methylsig(bsseq, group_column = group_column, 
    comparison_groups = c("case" = case, "control" = control),
    disp_groups = c("case" = TRUE, "control" = TRUE))
  
  # Add diff_meth_results as metadata to genomic_regions and return
  overlaps = data.frame(GenomicRanges::findOverlaps(genomic_regions, diff_meth_results, type = "equal"))
  S4Vectors::mcols(genomic_regions) = data.frame(matrix(NA, 
    nrow = length(genomic_regions), ncol = length(S4Vectors::mcols(diff_meth_results)), 
    dimnames = list(names(genomic_regions), names(S4Vectors::mcols(diff_meth_results)))))
  S4Vectors::mcols(genomic_regions)[overlaps$queryHits, ] = data.frame(S4Vectors::mcols(diff_meth_results)[overlaps$subjectHits, ])
  return(genomic_regions)
  
}
