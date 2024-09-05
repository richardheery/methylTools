#' Summarize methylation values for CpG sites in regions within a chunk
#' 
#' @param chunk_regions Chunk with genomic regions of interest. 
#' @param meth_rse A RangedSummarizedExperiment with methylation values.
#' @param assay The assay from meth_rse to extract values from. 
#' Should be either an index or the name of an assay.
#' @param summary_function_list A named list of functions to perform. 
#' @param na.rm TRUE or FALSE indicating whether to remove NA values when calculating summaries.
#' @return A list of data.frames with summaries of the methylation of CpG sites overlapping regions in the supplied chunk.
.summarize_cpgs_in_chunk = function(chunk_regions, meth_rse, genomic_regions, assay, summary_function_list, na.rm){
  
  # Subset meth_rse_for_chunk for regions overlapping chunk_regions
  meth_rse_for_chunk <- subsetByOverlaps(meth_rse, chunk_regions)
  invisible(gc()) 
  
  # Find the overlaps of chunk_regions and meth_rse_for_chunk
  overlaps_df <- data.frame(findOverlaps(chunk_regions, meth_rse_for_chunk))
  
  # Add region names and CpG names to overlaps_df
  overlaps_df$genomic_region_name <- names(chunk_regions)[overlaps_df$queryHits]
  
  # Create a list matching region names to rows of meth_rse_for_chunk
  region_names_to_rows_list <- split(overlaps_df$subjectHits, overlaps_df$genomic_region_name)
  
  # Read all values from specified assay of meth_rse_for_chunk into memory and run the garbage collection
  meth_values <- as.matrix(SummarizedExperiment::assay(meth_rse_for_chunk, i = assay))
  row.names(meth_values) <- as.character(rowRanges(meth_rse_for_chunk))
  gc()
  
  # Create a combined summary function
  combined_summary_function = function(data){
    df = data.frame(lapply(summary_function_list, function(sf)
      sf(data, na.rm = na.rm)))
    df = tibble::rownames_to_column(df, "cpg_name")
  }
  
  # Summarize methylation values 
  meth_summary <- lapply(region_names_to_rows_list, function(x) 
    combined_summary_function(meth_values[x, , drop = FALSE]))
  rm(meth_values); gc()
  return(meth_summary)
  
}

#' Summarize methylation values for CpG sites in regions
#' 
#' @param meth_rse A RangedSummarizedExperiment with methylation values.
#' @param genomic_regions GRanges object with regions to summarize methylation values for. 
#' @param assay The assay from meth_rse to extract values from. 
#' Should be either an index or the name of an assay. Default is the first assay. 
#' @param summary_function_list A named list of functions to perform. 
#' Should be one of the row or column summary functions from matrixStats
#' @param na.rm TRUE or FALSE indicating whether to remove NA values when calculating summaries. Default is TRUE.
#' @param BPPARAM A BiocParallelParam object. Defaults to `BiocParallel::bpparam()`.
#' @return A list of data.frames with summaries of the methylation of CpG sites overlapping each region in genomic_regions. 
summarize_cpgs_in_regions = function(meth_rse, genomic_regions, assay = 1, summary_function_list, na.rm = TRUE, BPPARAM = BiocParallel::bpparam()){
  
  # Check that summary_function_list is a named list of functions
  if(is.null(names(summary_function_list)) || !any(sapply(summary_function_list, function(x) is(x, "function")))){
    stop("summary_function_list should be a named list of functions")
  }
  
  # Split genomic regions into chunks based on the number of methylation sites that they cover
  genomic_region_bins <- .chunk_regions(meth_rse = meth_rse, genomic_regions = genomic_regions, ncores = BiocParallel::bpnworkers(BPPARAM))
  
  # For each sequence get methylation of the associated regions
  message("Summarizing genomic region methylation")
  cpg_summaries <- BiocParallel::bpmapply(FUN = .summarize_cpgs_in_chunk, 
    chunk_regions = genomic_region_bins, MoreArgs = list(meth_rse = meth_rse, 
      assay = assay, summary_function_list = summary_function_list, na.rm = na.rm), 
    BPPARAM = BPPARAM, SIMPLIFY = FALSE)
  
  # Flatten cpg_summaries, combine results into a single table and return 
  cpg_summaries = do.call(c, unname(cpg_summaries))
  cpg_summaries <- dplyr::bind_rows(cpg_summaries, .id = "genomic_region")
  return(cpg_summaries)
  
}