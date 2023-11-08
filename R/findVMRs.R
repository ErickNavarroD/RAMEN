#'  Map revmap column to probe names after reducing a GenomicRanges object
#'
#'  Given a revmap row (e.g. 1 5 6), we map those positions to their corresponding probe names
#'  (and end up with something like "cg00000029",  "cg00000158", "cg00000165".This is a helper function
#'  of findVMRs()).
#'
#' @param positions A revmap row in the form of a vector
#' @param manifest_hvp the manifest of the highly variable probes used in the findVMRs() function
#' with the probes as row names
#'
#' @return a vector with the names of the probes that conform one reduced region
#'
map_revmap_names = function(positions, manifest_hvp){
  #We start with 1 5 6
  #We want to end with cg00000029,  cg00000158 cg00000165
  names = c()
  for (element in positions){
    names =c(names, row.names(manifest_hvp)[element] )
  }
  return(names)
}



#' Identify Variable Methylated Regions in microarrays
#'
#' Identifies autosomal Highly Variable Probes (HVP) and merges them into Variable Methylated Regions (VMRs) given an Illumina manifest.
#'
#' This function identifies HVPs using MAD scores or variance metrics, and groups them into VMRs, which are defined as clusters of proximal and correlated HVPs (distance and correlation defined by the user). Output VMRs can be separated into canonical and non canonical. Canonical VMRs are regions that meet the correlation and closeness criteria. For guidance on which correlation threshold to use, we recommend checking the Supplementary Figure 1 of the CoMeBack R package (Gatev *et al.*, 2020) where a simulation to empirically determine a default guidance specification for a correlation threshold parameter dependent on sample size is done. As default, we use a threshold of 0.15 as per the CoMeBack authors minimum threshold suggestion. On the other hand, non canonical VMRs are regions that are composed of HVPs that have no nearby probes measured in the array (according to the max_distance parameter); this category was created to account for the Illumina EPIC array design, which has a high number of probes in regulatory regions that are represented by a single probe. Furthermore, these probes have been shown to be good representatives of the methylation state of its surroundings (Pidsley et al., 2016). By creating this category, we recover those informative HVPs that otherwise would be excluded from the analysis because of the array design.
#'
#' This function uses GenomicRanges::reduce() to group the regions, which is strand-sensitive. In the Illumina microarrays, the MAPINFO for all the probes
#' is usually provided as for the + strand. If you are using this array, we recommend to first
#' convert the strand of all the probes to "+".
#'
#' This function supports parallel computing for increased speed. To do so, you have to set the parallel backend
#' in your R session BEFORE running the function (e.g., doFuture::registerDoFuture()) and then the evaluation strategy (e.g., future::plan(multisession)). After that,
#' the function can be run as usual. When working with big datasets, the parallel backend might throw an error if you exceed
#' the maximum allowed size of globals exported for future expression. This can be fixed by increasing the allowed size (e.g. running options(future.globals.maxSize= +Inf) )
#'
#'Note: this function excludes sex chromosomes.
#'
#' @param array_manifest Information about the probes on the array. Requires the columns MAPINFO (basepair position
#' of the probe in the genome), CHR (chromosome), TargetID (probe name) and STRAND (this is very important to set up, since
#' the VMRs will only be created based on CpGs on the same strand; if the positions are reported based on a single DNA strand,
#' this should contain either a vector of only "+", "-" or "*" for all of the probes).
#' @param methylation_data A data frame containing M or B values, with samples as columns and probes as rows. Data is expected to have already passed through quality control and cleaning steps.
#' @param cor_threshold Numeric value (0-1) to be used as the median pearson correlation threshold for identifying VMRs (i.e.
#' all VMRs will have a median pairwise probe correlation of this parameter).
#' @param var_method Method to use to measure variability in the data set. The options are "mad" (median absolute deviation)
#' or "variance".
#' @param var_threshold_percentile The percentile (0-1) to be used as cutoff to define Highly Variable Probes (and
#' therefore VMRs). The default is 0.9 because this percentile has been traditionally used in previous studies.
#' @param max_distance Maximum distance allowed for two probes to be grouped into a region. The default is 1000
#' because this window has been traditionally used in previous studies.
#'
#' @return A list with the following elements:
#'  - $var_score_threshold: threshold used to define Highly Variable Probes (mad or variance, depending on the specified choice).
#'  - $highly_variable_probes: a data frame with the probes that passed the variability score threshold imposed by the user, and their variability score (MAD score or variance).
#'  - $canonical_VMRs: a GRanges object with strict candidate VMRs - regions composed of two or more
#'   contiguous, correlated and proximal Highly Variable Probes; thresholds depend on the ones specified
#'    by the user)
#'  - $non_canonical_VMRs: a GRanges object with highly variable probes without neighboring
#'  CpGs measured in *max_distance* on the array. Category created to take into acccount the Illumina array design of single probes capturing the methylation state of regulatory regions.
#'
#' @export
#'@examples
#' #We need to modify the RAMEN::test_array_manifest object by assigning to
#' #row names to the probe ID column; it was saved this way because storing
#' #the TargetID as row names reduced significantly the size of the data set.
#' test_array_manifest_final = RAMEN::test_array_manifest %>%
#' tibble::rownames_to_column(var = "TargetID")
#'
#' VMRs = RAMEN::findVMRs(array_manifest = test_array_manifest_final,
#'                        methylation_data = RAMEN::test_methylation_data,
#'                        cor_threshold = 0,
#'                        var_method = "variance",
#'                        var_threshold_percentile = 0.9,
#'                        max_distance = 1000)
#'
findVMRs = function(array_manifest,
                    methylation_data,
                    cor_threshold = 0.15,
                    var_method = "variance",
                    var_threshold_percentile = 0.9,
                    max_distance = 1000){
  #Check that the array manifest is in the right format
  if(!all(c("MAPINFO","CHR","TargetID","STRAND") %in% colnames(array_manifest))){
    stop("Please make sure the array manifest has the required columns with the appropiate names (check documentation for further information)")
  }
  #Check that the array strand is in the format expected by the user
  if(base::length(base::unique(array_manifest$STRAND)) > 1){
    warning("The manifest currently has more than one type of strands. Please note that this function is strand sensitive. So, probes in proximal coordinates but different strands on the manifest will not be grouped together. We recommend setting all of the probes to the same strand for arrays such as EPIC")
  }
  #Check that the method choice is correct
  if(var_method == "mad"){
    var_scores = apply(methylation_data, 1, stats::mad) %>%
      as.data.frame() %>%
      dplyr::rename("var_score" = ".")
  } else if (var_method == "variance") {
    var_scores = apply(methylation_data, 1, stats::var) %>%
      as.data.frame() %>%
      dplyr::rename("var_score" = ".")
  } else {
    stop("The method must be either 'mad' or 'variance'. Please select one of those options")
  }

  ####Identify highly variable probes ####
  message("Identifying Highly Variable Probes...")
  var_threshold = stats::quantile(var_scores$var_score, var_threshold_percentile)
  #Filter the manifest to remove the probes that have no variability score information because they were not measured/did not pass the QC and are not highly variable
  manifest_hvp = array_manifest %>%
    dplyr::select(c(TargetID, CHR, MAPINFO, STRAND)) %>%
    dplyr::filter(!is.na(MAPINFO), #Remove probes with no map info
                  !CHR %in% c("X","Y"), #Remove sexual chromosomes
                  TargetID %in% row.names(var_scores %>%
                                            dplyr::filter(var_score >= var_threshold))) %>% #Remove probes that have no methylation information in the processed data and are not highly variable
    dplyr::left_join(var_scores %>% #Add variability information
                       tibble::rownames_to_column(var = "TargetID"),
                     by = "TargetID") %>%
    dplyr::mutate(CHR = droplevels(CHR)) %>%
    dplyr::arrange(CHR) #important step for using Rle later when constructing the GenomicRanges object!
  rownames(manifest_hvp) = manifest_hvp$TargetID

  #### Identify probes with no neighbours####
  message("Identifying non canonical Variable Methylated Regions...")
  full_manifest = array_manifest %>%
    dplyr::select(c(TargetID, CHR, MAPINFO, STRAND)) %>%
    dplyr::filter(!is.na(MAPINFO), #Remove probes with no map info
                  !CHR %in% c("X","Y"), #Remove sexual chromosomes
                  TargetID %in% row.names(var_scores)) %>%  #keep only the probes where we have methylation information
    dplyr::mutate(CHR = droplevels(CHR)) %>%
    dplyr::arrange(CHR) #important step for using Rle later when constructing the GenomicRanges object!
  rownames(full_manifest) = full_manifest$TargetID

  #Convert the full manifest to a GenomicRanges object
  seqnames_full_manifest_gr = table(full_manifest$CHR)
  full_manifest_gr = GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(names(seqnames_full_manifest_gr), as.numeric(seqnames_full_manifest_gr)), #Number of chromosome; as.numeric to convert from table to numeric vector
    ranges = IRanges::IRanges(full_manifest$MAPINFO, end = full_manifest$MAPINFO ,
                              names = full_manifest$TargetID),
    strand = S4Vectors::Rle(rle(as.character(full_manifest$STRAND))$values,
                            rle(as.character(full_manifest$STRAND))$lengths ))

  #### Group the probes into regions to detect non-canonical VMRs
  regions_full_manifest = GenomicRanges::reduce(full_manifest_gr, with.revmap = TRUE, min.gapwidth = max_distance)
  #Add the number of probes in each region
  S4Vectors::mcols(regions_full_manifest)$n_probes = sapply(S4Vectors::mcols(regions_full_manifest)$revmap, length)
  #Substitute revmap with the name of the probes in each region
  S4Vectors::mcols(regions_full_manifest)$probes = sapply(S4Vectors::mcols(regions_full_manifest)$revmap, map_revmap_names, full_manifest)
  #Remove revmap mcol
  S4Vectors::mcols(regions_full_manifest)$revmap = NULL
  #Keep elements with only one probe
  lonely_probes = regions_full_manifest[(GenomicRanges::elementMetadata(regions_full_manifest)[,"n_probes"] <= 1)] %>%
    as.data.frame() %>%
    dplyr::pull(probes) %>%
    unlist()

  #### Identify VMRs####
  message("Identifying canonical Variable Methylated Regions...")
  #convert the highly variable probes data frame to a GenomicRanges object
  seqnames_gr = table(manifest_hvp$CHR)
  gr = GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(names(seqnames_gr), as.numeric(seqnames_gr)), #Number of chromosome; as.numeric to convert from table to numeric vector
    ranges = IRanges::IRanges(manifest_hvp$MAPINFO, end = manifest_hvp$MAPINFO ,
                              names = manifest_hvp$TargetID),
    strand = S4Vectors::Rle(rle(as.character(manifest_hvp$STRAND))$values,
                            rle(as.character(manifest_hvp$STRAND))$lengths ),
    var_score = manifest_hvp$var_score) #Metadata

  #Create the regions
  candidate_VMRs = GenomicRanges::reduce(gr, with.revmap = TRUE, min.gapwidth = max_distance)
  #Add the number of probes in each region
  S4Vectors::mcols(candidate_VMRs)$n_VMPs = sapply(S4Vectors::mcols(candidate_VMRs)$revmap, length)
  #Add the width of each region
  S4Vectors::mcols(candidate_VMRs)$width = S4Vectors::width(candidate_VMRs)
  #Substitute revmap with the name of the probes in each VMR
  S4Vectors::mcols(candidate_VMRs)$probes = sapply(S4Vectors::mcols(candidate_VMRs)$revmap, map_revmap_names, manifest_hvp)
  #Remove revmap mcol
  S4Vectors::mcols(candidate_VMRs)$revmap = NULL

  ### Capture canonical VMRs ###
  message("Applying correlation filter to canonical Variable Methylated Regions...")
  canonical_VMRs = candidate_VMRs[(GenomicRanges::elementMetadata(candidate_VMRs)[,"n_VMPs"] > 1)] %>%
    #Check for correlation between probes in these strict regions #
    as.data.frame() %>% #Convert the GR to a data frame so that I can use medCorVMR() and neigbouring_check()
    ### Check that the VMRs contain surrounding probes
    medCorVMR(VMR_df = ., methylation_data = methylation_data) %>% # Compute the median correlation of each region
    dplyr::filter(median_correlation > cor_threshold) %>%  #Remove VMRs whose CpGs are not correlated
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) #Create a GR object again
  colnames(S4Vectors::mcols(canonical_VMRs))[2] = "width" #Changing the name of one metadata variable that was modified when transforming from data frame to GR object

  ### Capture non-canonical VMRs ###
  non_canonical_VMRs =  candidate_VMRs[(GenomicRanges::elementMetadata(candidate_VMRs)[,"n_VMPs"] <= 1)] #Select the VMRs with 1 probe per region
  non_canonical_VMRs =  candidate_VMRs[(GenomicRanges::elementMetadata(candidate_VMRs)[,"probes"] %in% lonely_probes)] #Select the lonely probes
  GenomicRanges::mcols(non_canonical_VMRs)$median_correlation = rep(NA, nrow(GenomicRanges::mcols(non_canonical_VMRs))) #Add a column of NAs under the name of median_correlation to match the strict_VMRs



  return(list(
    var_score_threshold = var_threshold,
    highly_variable_probes = var_scores %>%
      tibble::rownames_to_column(var = "TargetID") %>%
      dplyr::filter(TargetID %in% manifest_hvp$TargetID),
    canonical_VMRs = canonical_VMRs,
    non_canonical_VMRs = non_canonical_VMRs
  ))
}

