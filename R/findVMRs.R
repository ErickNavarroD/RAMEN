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
#' Identifies Highly Variable Probes (HVP) and merges them into Variable Methylated Regions (VMRs) given an Illumina manifest.
#'
#' This function identifies HVPs using MAD scores or variance metrics, and groups them into VMRs, which are defined as clusters of proximal and correlated HVPs (distance and correlation defined by the user). To identify VMR, RAMEN::findVMRs() relies first on the identification of Highly Variable Probes in a data set. We support two methods for labelling probes as highly variable in the data set: 1)
#'
#'  Output VMRs can be separated into canonical and non canonical. Canonical VMRs are regions that meet the correlation and closeness criteria. For guidance on which correlation threshold to use, we recommend checking the Supplementary Figure 1 of the CoMeBack R package (Gatev *et al.*, 2020) where a simulation to empirically determine a default guidance specification for a correlation threshold parameter dependent on sample size is done. As default, we use a threshold of 0.15 as per the CoMeBack authors minimum threshold suggestion. On the other hand, non canonical VMRs are regions that are composed of HVPs that have no nearby probes measured in the array (according to the max_distance parameter); this category was created to account for the Illumina EPIC array design, which has a high number of probes in regulatory regions that are represented by a single probe. Furthermore, these probes have been shown to be good representatives of the methylation state of its surroundings (Pidsley et al., 2016). By creating this category, we recover those informative HVPs that otherwise would be excluded from the analysis because of the array design.
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
#'Note: this function does not exclude sex chromosomes. If you want to exclude them, you can do so in the methylation_data object before running the function.
#'
#' @param array_manifest Information about the probes on the array in a format compatible with the Bioconductor annotation packages. The user can specify one of the supported human microarrays ("IlluminaHumanMethylation450k" with the hg19 genome build, "IlluminaHumanMethylationEPICv1" with the hg19 genome build, or "IlluminaHumanMethylationEPICv2" with the hg38 genome build), or provide a manifest. The manifest requires the probe names as row names, and the following columns: "chr" (chromosome); "pos" (genomic location of the probe in the genome); and "strand" (this is very important to set up, since the VMRs will only be created based on CpGs on the same strand; if the positions are reported based on a single DNA strand, this should contain either a vector of only "+", "-" or "*" for all of the probes).
#' @param methylation_data A data frame containing M or B values, with samples as columns and probes as rows. Data is expected to have already passed through quality control and cleaning steps.
#' @param cor_threshold Numeric value (0-1) to be used as the median pearson correlation threshold for identifying VMRs (i.e.
#' all VMRs will have a median pairwise probe correlation higher than this threshold).
#' @param var_method A string indicating the method to use to measure variability in the data set. The options are "mad" (median absolute deviation)
#' or "variance".
#' @param var_distribution A string indicating which probes in the data set should be used to create the variability distribution, from which the variability threshold is taken from (percentile threshold determined by var_threshold_percentile). The options are "ultrastable" (a subset of CpGs that are stably methylated/unmethylated across human tissues and developmental states described by [Edgar R., et al.](https://doi.org/10.1186/1756-8935-7-28) in 2014); and "all" (all probes in the data set). The "ultrastable" option is only compatible with Illumina human microarrays. The default is "ultrastable".
#' @param var_threshold_percentile The percentile (0-1) to be used as cutoff to define Highly Variable Probes (which are then grouped into VMRs). If using the variability of the "ultrastable" probes, we recommend a high threshold (default is 0.99), since these probes are expected to display a very low variation in human tissues. If using the variability of "all" probes, we recommend using a percentile of 0.9 since it captures the top 10% most variable probes and has been traditionally used in previous studies.
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
#' @examples
#'
#' VMRs = RAMEN::findVMRs(methylation_data = RAMEN::test_methylation_data,
#'                        array_manifest = "IlluminaHumanMethylationEPICv1",
#'                        cor_threshold = 0.15,
#'                        var_method = "variance",
#'                        var_distribution = "ultrastable",
#'                        var_threshold_percentile = 0.99,
#'                        max_distance = 1000)
#'
findVMRs = function(methylation_data,
                    array_manifest,
                    cor_threshold = 0.15,
                    var_method = "variance",
                    var_distribution = "ultrastable",
                    var_threshold_percentile = 0.99,
                    max_distance = 1000){
  #Check that the array manifest is in the right format
  if(is.data.frame(array_manifest)){
    if(!all(c("chr","pos", "strand") %in% colnames(array_manifest))) stop("The array_manifest data frame does not have the required columns. Please provide a manifest with the required columns or provide a string with one of the supported human microarrays ('IlluminaHumanMethylation450k', 'IlluminaHumanMethylationEPICv1','IlluminaHumanMethylationEPICv2')")
    #Check that the array strand is in the format expected by the user
    if(base::length(base::unique(array_manifest$strand)) > 1) warning("The manifest currently has more than one type of strands. Please note that this function is strand sensitive. So, probes in proximal coordinates but different strands on the manifest will not be grouped together. Many array manifests such as the Illumina EPIC one include the PROBE strand, but the position of the actual CpGs (pos) is reported in the same strand; in those cases we recommend setting all of the probes to the same strand.")
    if(var_distribution == "ultrastable") {
      #If the user provides their own manifest and is choosing to use the ultrastable probes, make sure that a good number of them is present in the data set. If not, throw an error
      if(sum(row.names(array_manifest) %in% RAMEN::ultrastable_cpgs) < 100) stop ("The var_distribution = 'ultrastable' option is only compatible with Illumina human microarrays at the moment. If you are using a human Illumina microarray please indicate it with their corresponding string, or make sure that it contains a more than 100 ultrastable probes (RAMEN::ultrastable_cpgs). If not, please get the variability threshold based on all the probes in your data set(var_distribution = 'all', var_threshold_percentile = 0.9). ")
    } else if(is.character(array_manifest)){
    }
    if(!array_manifest %in% c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPICv1","IlluminaHumanMethylationEPICv2"))stop("The string you provided in array_manifest is not currently supported in RAMEN. Please provide a manifest with the required columns or provide a string with one of the supported human microarrays ('IlluminaHumanMethylation450k', 'IlluminaHumanMethylationEPICv1','IlluminaHumanMethylationEPICv2')")
  } else {
    stop("The array_manifest object is not a data.frame nor a string. Please provide a manifest with the required columns or provide a string with one of the supported human microarrays ('IlluminaHumanMethylation450k', 'IlluminaHumanMethylationEPICv1','IlluminaHumanMethylationEPICv2')")
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
  # Get the variability threshold
  if(var_distribution == "all") {
    var_threshold = stats::quantile(var_scores$var_score,
                                    var_threshold_percentile)
  } else if (var_distribution == "ultrastable"){
    if(array_manifest == "IlluminaHumanMethylationEPICv2"){
      #Get the name of the ultrastable probes in the EPICv2 format
      epicv2_ultrastable_cpgs = IlluminaHumanMethylationEPICv2anno.20a1.hg38::Other |>
        data.frame() |>
        dplyr::filter(Methyl450_Loci %in% RAMEN::ultrastable_cpgs) |>
        tibble::rownames_to_column("epicv2_probes") |>
        dplyr::pull(epicv2_probes)
      var_threshold = stats::quantile(var_scores[(row.names(var_scores) %in% epicv2_ultrastable_cpgs),],
                                      var_threshold_percentile)
    } else {
      #EPICv1 or 450k (same probe name as the ultrastable probes)
      var_threshold = stats::quantile(var_scores[(row.names(var_scores) %in% RAMEN::ultrastable_cpgs),],#Subset only ultrastable probes
                                      var_threshold_percentile)
    }
  }
  #Replace the array manifest if the user provided a string with the name of the array
  if(is.character(array_manifest)){
    if(array_manifest == "IlluminaHumanMethylation450k"){
      manifest = data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations)
    } else if(array_manifest == "IlluminaHumanMethylationEPICv1"){
      manifest = data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations)
    } else if(array_manifest == "IlluminaHumanMethylationEPICv2"){
      manifest = data.frame(IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations)
    }
  } else manifest = array_manifest
  #Filter the manifest to remove the probes that have no variability score information because they were not measured/did not pass the QC and are not highly variable
  manifest_hvp = manifest %>%
    tibble::rownames_to_column(var = "TargetID") %>%
    dplyr::select(c(TargetID, chr, pos, strand)) %>%
    dplyr::filter(!is.na(pos), #Remove probes with no map info
                  TargetID %in% row.names(var_scores %>%
                                            dplyr::filter(var_score >= var_threshold))) %>% #Remove probes that have no methylation information in the processed data and are not highly variable
    dplyr::left_join(var_scores %>% #Add variability information
                       tibble::rownames_to_column(var = "TargetID"),
                     by = "TargetID") %>%
    dplyr::arrange(chr) %>%  #important step for using Rle later when constructing the GenomicRanges object!
    as.data.frame()
  rownames(manifest_hvp) = manifest_hvp$TargetID
  if(is.factor(manifest_hvp$chr)) manifest_hvp = manifest_hvp %>% dplyr::mutate(chr = droplevels(chr))

  #### Identify probes with no neighbours####
  message("Identifying non canonical Variable Methylated Regions...")
  full_manifest = manifest %>%
    tibble::rownames_to_column(var = "TargetID") %>%
    dplyr::select(c(TargetID, chr, pos, strand)) %>%
    dplyr::filter(!is.na(pos), #Remove probes with no map info
                  TargetID %in% row.names(var_scores)) %>%  #keep only the probes where we have methylation information
    dplyr::arrange(chr) %>% #important step for using Rle later when constructing the GenomicRanges object!
    as.data.frame()
  rownames(full_manifest) = full_manifest$TargetID
  if(is.factor(full_manifest$chr)) full_manifest = full_manifest %>% dplyr::mutate(chr = droplevels(chr))

  #Convert the full manifest to a GenomicRanges object
  seqnames_full_manifest_gr = table(full_manifest$chr)
  full_manifest_gr = GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(names(seqnames_full_manifest_gr), as.numeric(seqnames_full_manifest_gr)), #Number of chromosome; as.numeric to convert from table to numeric vector
    ranges = IRanges::IRanges(full_manifest$pos, end = full_manifest$pos ,
                              names = full_manifest$TargetID),
    strand = S4Vectors::Rle(rle(as.character(full_manifest$strand))$values,
                            rle(as.character(full_manifest$strand))$lengths ))

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
  seqnames_gr = table(manifest_hvp$chr)
  gr = GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(names(seqnames_gr), as.numeric(seqnames_gr)), #Number of chromosome; as.numeric to convert from table to numeric vector
    ranges = IRanges::IRanges(manifest_hvp$pos, end = manifest_hvp$pos ,
                              names = manifest_hvp$TargetID),
    strand = S4Vectors::Rle(rle(as.character(manifest_hvp$strand))$values,
                            rle(as.character(manifest_hvp$strand))$lengths ),
    var_score = manifest_hvp$var_score) #Metadata

  #Create the regions
  candidate_VMRs = GenomicRanges::reduce(gr, with.revmap = TRUE, min.gapwidth = max_distance)
  #Add the number of probes in each region
  S4Vectors::mcols(candidate_VMRs)$n_VMPs = sapply(S4Vectors::mcols(candidate_VMRs)$revmap, length)
  #Add the width of each region
  #S4Vectors::mcols(candidate_VMRs)$width = S4Vectors::width(candidate_VMRs)
  #Substitute revmap with the name of the probes in each VMR
  S4Vectors::mcols(candidate_VMRs)$probes = sapply(S4Vectors::mcols(candidate_VMRs)$revmap, map_revmap_names, manifest_hvp)
  #Remove revmap mcol
  S4Vectors::mcols(candidate_VMRs)$revmap = NULL

  ### Capture canonical VMRs ###
  message("Applying correlation filter to canonical Variable Methylated Regions...")
  canonical_VMRs = candidate_VMRs[(GenomicRanges::elementMetadata(candidate_VMRs)[,"n_VMPs"] > 1)] %>%
    #Check for correlation between probes in these strict regions #
    as.data.frame() #Convert the GR to a data frame so that I can use medCorVMR() and neigbouring_check()
  ### Check that the VMRs contain surrounding probes only if we have potential canonical VMRs
  if(nrow(canonical_VMRs) > 0){
    canonical_VMRs = canonical_VMRs %>%
      medCorVMR(VMR_df = ., methylation_data = methylation_data) %>% # Compute the median correlation of each region
      dplyr::filter(median_correlation > cor_threshold) %>%  #Remove VMRs whose CpGs are not correlated
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) #Create a GR object again
  } else warning("No canonical VMRs were found in this data set")

  ### Capture non-canonical VMRs ###
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


