#' Convert Bismark genome-wide cytosine reports to BedGraph files with methylation values for each CpG site.
#' 
#' @param cytosine_report_files Paths to genome-wide cytosine reports output by bismark_methylation_extractor --cytosine_report.
#' @param files_0_based A logical value indicating if files are 0-based. Default is FALSE.
#' @param genome A BSgenome object or the name of a BSgenome package from which to extract CpG sites.
#' @param require_both_strands A logical value indicating whether to remove CpG sites in a file
#' without reads coming from both strands. Default is TRUE. 
#' @param output_directory Directory to save output bedGraphs. Directory is created if it doesn't already exist.
#' @param n_parallel_files Number of files to process in parallel. Default is 1.
#' @param decimal_places Number of decimal places to round methylation values to. Default value is 2.  
#' @return A vector with the paths to the created bedGraph files.
#' @export
convert_bismarck_cytosine_reports_to_bedgraphs = function(cytosine_report_files, files_0_based = FALSE, 
  genome, require_both_strands = TRUE, output_directory, n_parallel_files = 1, decimal_places = 2){
  
  # Extract CpG sites from indicated genome
  genome = BSgenome::getBSgenome(genome)
  message(paste("Extracting CpG sites from", attributes(genome)$pkgname))
  cpg_sites = extract_cpgs_from_bsgenome(genome)
  message(paste("Extracted", length(cpg_sites), "CpG sites from", attributes(genome)$pkgname))
    
  # Get seqlevelsStyle of cpg_sites from genome
  cpg_sites_seqlevels_style = seqlevelsStyle(getBSgenome(genome))
  
  # Create separate GRanges objects with the C and G position of CpG sites
  c_sites = resize(cpg_sites, width = 1, fix = "start")
  g_sites = IRanges::shift(c_sites, 1)
  
  # Create output directory if it doesn't already exist.
  if(!dir.exists(output_directory)){
    dir.create(output_directory)
  }
    
  # Get the basenames of the cytosine_report_files
  basenames = basename(gsub("\\..*", "", cytosine_report_files))
  
  # Create paths for output file
  output_files = paste0(output_directory, "/", basenames, ".bg.gz")
  
  # Create a cluster and register it if n_parallel_files is greater than 1
  if(n_parallel_files > 1){
    cluster = parallel::makeCluster(n_parallel_files, outfile = "")
    doParallel::registerDoParallel(cluster)
    `%dopar%` = foreach::`%dopar%`
    on.exit(parallel::stopCluster(cluster))
  } else {
    `%dopar%` = foreach::`%do%`
  }
  
  # Loop through each input file and process it
  foreach::foreach(in_file = cytosine_report_files, out_file = output_files, .packages = "GenomicRanges") %dopar% {
    
    # Print the name of the file being processed
    message(paste("Starting processing", in_file))
    
    # Read in the input file only selecting the first five columns and name them
    cytosine_report = data.table::fread(in_file, nThread = 1, select = 1:5, 
      col.names = c("seqnames", "start", "strand", "meth_reads", "unmeth_reads"), showProgress = F)
    
    # Convert the table into a GRanges object
    cytosine_report = GenomicRanges::makeGRangesFromDataFrame(cytosine_report, start.field = "start", 
      end.field = "start", keep.extra.columns = T)
    
    # Set seqlevelsStyle to that of cpg_sites
    GenomeInfoDb::seqlevelsStyle(cytosine_report) = cpg_sites_seqlevels_style
    
    # Create copies of c_sites and g_sites sites to store the number of 
    # methylated and unmethylated reads from the input file
    c_sites_temp = c_sites
    c_sites_temp$meth_value = NA
    c_sites_temp$unmeth_value = NA
    g_sites_temp = g_sites
    g_sites_temp$meth_value = NA
    g_sites_temp$unmeth_value = NA
    
    # Find overlaps between loci in the input file and C and G positions of CpG sites
    c_overlaps = data.frame(GenomicRanges::findOverlaps(c_sites_temp, cytosine_report))
    g_overlaps = data.frame(GenomicRanges::findOverlaps(g_sites_temp, cytosine_report))
    
    # Add the number of methylated and unmethylated reads to C and G positions of CpG sites from input file
    c_sites_temp$meth_value[c_overlaps$queryHits] = cytosine_report$meth_reads[c_overlaps$subjectHits]
    c_sites_temp$unmeth_value[c_overlaps$queryHits] = cytosine_report$unmeth_reads[c_overlaps$subjectHits]
    g_sites_temp$meth_value[g_overlaps$queryHits] = cytosine_report$meth_reads[g_overlaps$subjectHits]
    g_sites_temp$unmeth_value[g_overlaps$queryHits] = cytosine_report$unmeth_reads[g_overlaps$subjectHits]
    
    # Calculate the beta value for each CpG site as the sum of methylated reads for 
    # both the C and G position dividied by the total number of reads.
    # If require_both_strands is TRUE, will return NA if counts are missing for one of the strands
    c_sites_temp$beta_value = rowSums(cbind(c_sites_temp$meth_value, g_sites_temp$meth_value), na.rm = !require_both_strands)/
      rowSums(cbind(c_sites_temp$meth_value, g_sites_temp$meth_value, 
        c_sites_temp$unmeth_value, g_sites_temp$unmeth_value), na.rm = !require_both_strands)
    
    # Convert c_sites_temp into a data.frame, selecting only the position column and beta_value
    results_df = data.frame(c_sites_temp)[c("seqnames", "start", "end", "beta_value")]
    
    # Subtract 1 from the start so that it will be 0-based
    results_df$start = results_df$start - 1
    
    # Remove CpG sites where beta values are missing and round beta values to specified number of decimal places
    results_df = dplyr::filter(results_df, !is.na(beta_value))
    results_df$beta_value = round(results_df$beta_value, decimal_places)
    
    # Write results_df to indicated 
    data.table::fwrite(results_df, out_file, col.names = F, row.names = F, 
      quote = F, sep = "\t", nThread = 1, showProgress = F)
    
    # Return the number of CpGs for which beta values could be calculated
    message(paste("Finished processing", paste0(in_file, ":"), "Methylation levels calculated for", nrow(results_df), "CpG sites."))
    
  }
  
  # Return the paths to the output files
  return(output_files)
  
}

#' Create a GRanges with CpG sites from a BSgenome. 
#'
#' @param genome A BSgenome object or the name of a BSgenome package.
#' @param standard_sequences_only TRUE or FALSE indicating whether to only return sites 
#' on standard sequences (those without "-" in their names). Default is TRUE. 
#' @return A GRanges object CpG sites for the indicated genome. 
#' @export
#' @examples 
#' # Get human CpG sites for hg38 genome build
#' hg38_cpgs <- methodical::extractMethSitesFromGenome("BSgenome.Hsapiens.UCSC.hg38")
#' 
#' # Find CHG sites in Arabidopsis thaliana
#' arabidopsis_cphpgs <- methodical::extractMethSitesFromGenome("BSgenome.Athaliana.TAIR.TAIR9", pattern = "CHG")
extract_cpgs_from_bsgenome <- function(genome, standard_sequences_only = TRUE){
  
  # Check that inputs have the correct data type
  stopifnot(is(genome, "character") | is(genome, "BSgenome"), 
    S4Vectors::isTRUEorFALSE(standard_sequences_only))
  
  # If genome is a character, try to load genome with that name
  if(is.character(genome)){genome <- BSgenome::getBSgenome(genome)}
  
  # Find CpG sites in genome 
  cpg_gr <- GRanges(Biostrings::vmatchPattern("CG", genome, fixed = "subject"))
  
  # Filter for matches on "+" strand and set strand as "*"
  cpg_gr <- cpg_gr[GenomicRanges::strand(cpg_gr) == "+"]
  GenomicRanges::strand(cpg_gr) <- "*"
  
  # Add seqinfo to GRanges
  GenomeInfoDb::seqinfo(cpg_gr) <- GenomeInfoDb::seqinfo(genome)
  
  # Subset for standard sequences if specified
  if(standard_sequences_only){
    standard_sequences <- grep("_", names(genome), invert = TRUE, value = TRUE)
    GenomeInfoDb::seqlevels(cpg_gr, pruning.mode = "coarse") <- standard_sequences
  }
  
  # Resize cpg_gr so that it covers just the C position
  cpg_gr <- resize(cpg_gr, 1, fix = "start")
  
  # Sort and return cpg_gr
  return(sort(cpg_gr, ignore.strand = TRUE))
  
}

merge_bedgraph_c_and_g_values = function(input_bedgraphs, genome, output_directory, ncores){
  
  # Check that allowed values provided for genome
  match.arg(arg = genome, choices = c("hg19", "hg38"), several.ok = F)
  
  # Check that output directory doesn't already exist and create it otherwise
  if(dir.exists(output_directory)){
    stop("output_directory already exists")
  } else {
    dir.create(output_directory)
  }
  
  # Load CpG sites for appropriate genome
  if(genome == "hg19"){
    c_sites = methodicalFinal:::cpg_genome_ranges_hg19
  } else if(genome == "hg38"){
    c_sites = methodicalFinal:::cpg_genome_ranges_hg38
  }
  
  # Create GRanges with G sites of CpGs
  g_sites = IRanges::shift(c_sites, 1)
  
  `%dopar%` = foreach::`%dopar%`
  foreach::foreach(bg_in = input_bedgraphs) %dopar% {
    
    message(paste("Starting", bg_in))
    
    bg_name = basename(bg_in)
    bg_out = paste(output_directory, bg_name, sep = "/")
    bg_in = rtracklayer::import.bedGraph(bg_in)
  
    c_sites_temp = c_sites
    c_sites_temp$value = NA
    g_sites_temp = g_sites
    g_sites_temp$value = NA
    
    c_overlaps = data.frame(GenomicRanges::findOverlaps(c_sites_temp, bg_in))
    g_overlaps = data.frame(GenomicRanges::findOverlaps(g_sites_temp, bg_in))
    
    c_sites_temp$value[c_overlaps$queryHits] = bg_in$score[c_overlaps$subjectHits]
    g_sites_temp$value[g_overlaps$queryHits] = bg_in$score[g_overlaps$subjectHits]
    
    c_sites_temp$value = (c_sites_temp$value + g_sites_temp$value)/2
    results_df = data.frame(c_sites_temp)[c(1, 2, 3, 6)]
    results_df$start = results_df$start - 1
    results_df = dplyr::filter(results_df, !is.na(value))
    
    data.table::fwrite(results_df, bg_out, col.names = F, row.names = F, quote = F, sep = "\t", scipen = 4)
  }
  
}