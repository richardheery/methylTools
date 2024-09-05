#' Split genomic regions into balanced chunks based on the number of methylation sites that they cover
#' 
#' @param meth_rse A RangedSummarizedExperiment with methylation values.
#' @param genomic_regions A GRanges object.
#' @param max_sites_per_chunk The maximum number of methylation sites to load into memory at once for each chunk. 
#' @param ncores The number of cores that will be used. 
#' @return A GRangesList where each GRanges object overlaps approximately the number of methylation sites given by max_sites_per_chunk 
.chunk_regions <- function(meth_rse, genomic_regions, max_sites_per_chunk = floor(62500000/ncol(meth_rse)), ncores = 1){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_rse, "RangedSummarizedExperiment"), 
    is(genomic_regions, "GRanges"), 
    is(max_sites_per_chunk, "numeric") & max_sites_per_chunk >= 1 | is.null(max_sites_per_chunk),
    is(ncores, "numeric") & ncores >= 1)
  
  # Make seqlevels of genomic_regions the same as meth_rse
  GenomeInfoDb::seqlevels(genomic_regions, pruning.mode = "coarse") <- GenomeInfoDb::seqlevels(meth_rse)
  
  # Sort genomic_regions ordering by seqlevels of meth_rse
  genomic_regions <- sort(genomic_regions, ignore.strand = TRUE)
  
  # Get methylation sites covered in meth_rse
  meth_sites <- sort(rowRanges(meth_rse))
  
  # Find subset of meth_sites overlapping genomic_regions and the number of overlapping sites
  overlap_meth_sites <- subsetByOverlaps(meth_sites, genomic_regions)
  number_overlap_meth_sites <- length(overlap_meth_sites)
  
  # Find the number of chunks necessary
  nchunks <- ceiling(number_overlap_meth_sites/max_sites_per_chunk)
  
  # Ensure number of chunks is at least equal to the number of cores
  if(nchunks < ncores){nchunks <- ncores}
  
  # Update max_sites_per_chunk to be number_overlap_meth_sites divided by nchunks
  max_sites_per_chunk <- ceiling(number_overlap_meth_sites/nchunks)
  
  # Divide the overlapping methylation sites into bins
  overlap_meth_sites_bins <- ceiling(seq_along(overlap_meth_sites)/max_sites_per_chunk)
  
  # Create a GRangesList object with the methylation sites in each bin
  overlap_meth_sites_grl <- GenomicRanges::GRangesList(split(overlap_meth_sites, overlap_meth_sites_bins))
  
  # Find which regions overlap each bin
  genomic_region_bins <- lapply(seq_along(overlap_meth_sites_grl), function(x) 
    data.frame(bin = x, region = which(IRanges::overlapsAny(genomic_regions, overlap_meth_sites_grl[[x]]))))
  genomic_region_bins <- dplyr::bind_rows(genomic_region_bins)
  
  # Remove regions from a bin if they are already present in another bin. 
  genomic_region_bins <- genomic_region_bins[!duplicated(genomic_region_bins$region), ]
  
  # Turn region_bins into a list with the regions in each chunk and return
  genomic_region_bins <- GenomicRanges::GRangesList(split(genomic_regions[genomic_region_bins$region], 
    genomic_region_bins$bin))
  return(genomic_region_bins) 

}
