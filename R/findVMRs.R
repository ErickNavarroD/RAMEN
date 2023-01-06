#' Map revmap column to probe names after reducing a GenomicRanges object
#'
#'  Given a revmap row (e.g. 1 5 6), we map those positions to their corresponding probe names
#'  (and end up with something like "cg00000029",  "cg00000158", "cg00000165".This is a helper function
#'  of findVMRs().
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


#' Define VMRs using median absolute deviation (MAD) scores in microarrays
#'
#' Given a sorted illumina manifest and their MAD scores previously computed, identifies Variably
#' Methylated Probes (VMPs) and merge them into regions.
#' Note: this function uses GenomicRanges::reduce() to group the regions, which is strand-sensitive. In
#' the illumina microarrays, the MAPINFO for all the probes is usually provided as for the + strand. If
#' you are using this array, we recommend to previously convert the strand of all the probes to "+".
#'
#' @param array_manifest Information about the probes contained in the array. Must have the columns MAPINFO (position
#' of the probe in the genome), CHR (chromosome), TargetID(probe name) and STRAND (this is very important to set up, since
#' the VMRs will only be created on CpGs on the same strand; if the positions are reported based on a single DNA strand,
#' this should contain either a vector of only "+", "-" or "*" for all of the probes).
#' @param methylation_data A data frame containing M or B values, with samples as columns and probes as rows
#' @param MAD_threshold_percentile The MAD score percentile to be used as cut off to define Highly Variable Probes (and
#' therefore VMRs)
#' @param max_distance Maximum distance allowed for two probes to me grouped into a region
#' @param cor_threshold Numeric value (0-1) to be used as the median correlation threshold for identifying VMRs (i.e.
#' all VMRs will have a median pairwise probe correlation of this parameter).
#'
#' @return a list with the following elements:
#'  - $MAD_score_threshold: Threshold used to define highly variable probes.
#'  - $highly_variable_probes a data frame with the probes that passed the MAD score threshold imposed by the user, and their MAD score
#'  - $candidate_VMRs_strict: a GRanges object with strict candidate VMRs - regions composed of 2 or more
#'   highly variable probes closer thank 1 kb)
#'  - $candidate_VMRs_lonely_probes: a GRanges object with highly variable probes that had no neighboring
#'  CpGs measured in < 1kb in the array.
#'  @export
findVMRs = function(array_manifest, methylation_data, MAD_threshold_percentile = 0.9, max_distance = 1000,
                    cor_threshold){
  MAD_scores = apply(methylation_data, 1, stats::mad) %>%
    as.data.frame() %>%
    dplyr::rename("MAD_score" = ".")
  ####Identify highly variable probes ####
  MAD_threshold = stats::quantile(MAD_scores$MAD_score, MAD_threshold_percentile)
  #Filter the manifest to remove the probes that have no MAD score information because they were not measured/did not pass the QC and are not highly variable
  manifest_hvp = array_manifest %>%
    dplyr::select(c(TargetID, CHR, MAPINFO, STRAND)) %>%
    dplyr::filter(!is.na(MAPINFO), #Remove probes with no map info
                  !CHR %in% c("X","Y"), #Remove sexual chromosomes
                  TargetID %in% row.names(MAD_scores %>%
                                            dplyr::filter(MAD_scores >= MAD_threshold))) %>% #Remove probes that have no methylation information in the processed data and are not highly variable
    dplyr::left_join(MAD_scores %>% #Add MADscore information
                       tibble::rownames_to_column(var = "TargetID"), by = "TargetID") %>%
    dplyr::mutate(CHR = droplevels(CHR)) %>%
    dplyr::arrange(CHR)
  rownames(manifest_hvp) = manifest_hvp$TargetID

  #### Identify probes with no neighbours####
  full_manifest = array_manifest %>%
    dplyr::select(c(TargetID, CHR, MAPINFO, STRAND)) %>%
    dplyr::filter(!is.na(MAPINFO), #Remove probes with no map info
                  !CHR %in% c("X","Y"), #Remove sexual chromosomes
                  TargetID %in% row.names(MAD_scores)) %>%  #keep only the probes where we have methylation information
    dplyr::mutate(CHR = droplevels(CHR)) %>%
    dplyr::arrange(CHR)
  rownames(full_manifest) = full_manifest$TargetID

  #Convert the full manifest to a GenomicRanges object
  seqnames_full_manifest_gr = table(full_manifest$CHR)
  full_manifest_gr = GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(names(seqnames_full_manifest_gr), as.numeric(seqnames_full_manifest_gr)), #Number of chromosome
    ranges = IRanges::IRanges(full_manifest$MAPINFO, end = full_manifest$MAPINFO ,
                              names = full_manifest$TargetID),
    strand = S4Vectors::Rle(rle(as.character(full_manifest$STRAND))$values,
                            rle(as.character(full_manifest$STRAND))$lengths ))

  #Group the probes into regions
  regions_full_manifest = GenomicRanges::reduce(full_manifest_gr, with.revmap = TRUE, min.gapwidth = max_distance)
  #Add the number of probes in each region
  S4Vectors::mcols(regions_full_manifest)$n_probes = sapply(S4Vectors::mcols(regions_full_manifest)$revmap, length)
  #Substitute revmap with the name of the probes in each VMR
  S4Vectors::mcols(regions_full_manifest)$probes = sapply(S4Vectors::mcols(regions_full_manifest)$revmap, map_revmap_names, full_manifest)
  #Remove revmap mcol
  S4Vectors::mcols(regions_full_manifest)$revmap = NULL
  #Keep elements with only one probe
  lonely_probes = regions_full_manifest[(GenomicRanges::elementMetadata(regions_full_manifest)[,"n_probes"] <= 1)] %>%
    as.data.frame() %>%
    dplyr::pull(probes) %>%
    unlist()

  #### Identify candidate VMRs####
  #convert the highly variable probes data frame to a GenomicRanges object
  seqnames_gr = table(manifest_hvp$CHR)
  gr = GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(names(seqnames_gr), as.numeric(seqnames_gr)), #Number of chromosome
    ranges = IRanges::IRanges(manifest_hvp$MAPINFO, end = manifest_hvp$MAPINFO ,
                              names = manifest_hvp$TargetID),
    strand = S4Vectors::Rle(rle(as.character(manifest_hvp$STRAND))$values,
                            rle(as.character(manifest_hvp$STRAND))$lengths ),
    MAD_score = manifest_hvp$MAD_score) #Metadata

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

  ### Capture strict VMRs ###
  candidate_VMRs_strict = candidate_VMRs[(GenomicRanges::elementMetadata(candidate_VMRs)[,"n_VMPs"] > 1)] %>%
    #Check for correlation between probes in these strict regions #
    as.data.frame() %>% #Convert the GR to a data frame so that I can use medCorVMR() and neigbouring_check()
    ### Check that the VMRs contain surrounding probes
    medCorVMR(VMR_df = ., data_methylation = methylation_data) %>% # Compute the median correlation of each region
    #neighbFilterVMR() %>%
    dplyr::filter(median_correlation > cor_threshold) %>%  #Remove VMRs whose CpGs are not correlated
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) #Create a GR object again
  colnames(S4Vectors::mcols(candidate_VMRs_strict))[2] = "width" #Changing the name of one metadata variable that was modified when transforming from data frame to GR object

  ### Capture lonely probes VMRs ###
  candidate_VMRs_lonely_probes =  candidate_VMRs[(GenomicRanges::elementMetadata(candidate_VMRs)[,"n_VMPs"] <= 1)] #Select the VMRs with 1 probe per region
  candidate_VMRs_lonely_probes =  candidate_VMRs[(GenomicRanges::elementMetadata(candidate_VMRs)[,"probes"] %in% lonely_probes)] #Select the lonely probes
  GenomicRanges::mcols(candidate_VMRs_lonely_probes)$median_correlation = rep(NA, nrow(GenomicRanges::mcols(candidate_VMRs_lonely_probes))) #Add a column of NAs under the name of median_correlation to match the strict_VMRs



  return(list(
    MAD_score_threshold = MAD_threshold,
    highly_variable_probes = MAD_scores %>%
      tibble::rownames_to_column(var = "TargetID") %>%
      dplyr::filter(TargetID %in% manifest_hvp$TargetID),
    candidate_VMRs_strict = candidate_VMRs_strict,
    candidate_VMRs_lonely_probes = candidate_VMRs_lonely_probes
  ))
}

