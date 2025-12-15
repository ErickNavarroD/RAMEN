#'  Map revmap column to probe names after reducing a GenomicRanges object
#'
#'  Given a revmap row (e.g. 1 5 6), we map those positions to their corresponding probe names
#'  (and end up with something like "cg00000029",  "cg00000158", "cg00000165".This is a helper function of findVML()).
#'
#' @param positions A revmap row in the form of a vector
#' @param manifest_hvp the manifest of the highly variable probes used in the findVML() function
#' with the probes as row names
#'
#' @return a vector with the names of the probes that conform one reduced region
#'
#' @examples
#' \dontrun{
#'   target = data.frame(row.names = c("a", "b", "c", "d"), values = c(1,1,1,1))
#'   query = c(2,1)
#'
#'   map_revmap_names(positions = query, manifest_hvp = target)
#'   ## Expected output: c("b", "a")
#' }

map_revmap_names <- function(positions, manifest_hvp) {
  # We start with 1 5 6
  # We want to end with cg00000029,  cg00000158 cg00000165
  names <- c()
  for (element in positions) {
    names <- c(names, row.names(manifest_hvp)[element])
  }
  return(names)
}


#' Identify Variable Methylated Loci in microarrays
#'
#' Identifies Highly Variable Probes (HVP) and groups them into Variable Methylated Loci (VML) given an Illumina manifest.The output of this function provides the HVPs, and the identified VML, which are made of Variable Methylated Regions and sparse Variable Methylated Probes. See Details below for more information.
#'
#' This function identifies HVPs based on MAD scores or variance, and groups them into VML, which are defined as genomic regions with high DNA methylation variability.To best capture methylome variability patterns in microarrays, we identify two types of VML: Variably Methylated Regions (VMRs) and sparse Variably Methylated Probes (sVMPs) .
#'
#' In one hand, we defined VMRs as two or more proximal highly variable probes (default: < 1kb apart) with correlated DNAme level (default: r > 0.15). Modelling DNAme variability through regions rather than individual CpGs provides several methodological advantages in association studies, since CpGs display a significant correlation for co-methylation when they are close (less than or equal to 1 kilobase). Modelling DNAme variability through regions rather than individual CpGs provides several methodological advantages in association studies, since CpGs display a significant correlation for co-methylation when they are close (less than or equal to 1 kilobase)
#'
#' In addition to traditional VMRs, we also identified sparse Variably Methylated Probes (sVMPs), a second type of VML that takes into account the sparse and non-uniformly distributed coverage of CpGs in microarrays to tailor our analysis to this DNAme platform. sVMPs aimed to retain genomic regions with high DNAme variability measured by single probes, where probe grouping based on proximity and correlation is therefore not applicable. This is particularly relevant in the Illumina EPIC v1 array, where most covered regulatory regions (up to 93%) are represented by just one probe. Notably, based on empirical comparisons with whole-genome bisulfite sequencing data, these single probes are mostly representative of local regional DNAme levels due to their positioning (98.5-99.5%)
#'
#' This function uses GenomicRanges::reduce() to group the regions, which is strand-sensitive. In the Illumina microarrays, the MAPINFO for all the probes is usually provided for the + strand. If you are using this array, we recommend to first convert the strand of all the probes to "+".
#'
#' This function supports parallel computing for increased speed. To do so, you have to set the parallel backend
#' in your R session BEFORE running the function (e.g., *doParallel::registerDoParallel(4)*). After that, the function can be run as usual. When working with big datasets, the parallel backend might throw an error if you exceed the maximum allowed size of globals exported for future expression. This can be fixed by increasing the allowed size (e.g. running *options(future.globals.maxSize= +Inf)*)
#'
#' Note: this function does not exclude sex chromosomes. If you want to exclude them, you can do so in the methylation_data object before running the function.
#'
#' @param array_manifest Information about the probes on the array in a format compatible with the Bioconductor annotation packages. The user can specify one of the supported human microarrays ("IlluminaHumanMethylation450k" with the hg19 genome build, "IlluminaHumanMethylationEPICv1" with the hg19 genome build, or "IlluminaHumanMethylationEPICv2" with the hg38 genome build), or provide a manifest. The manifest requires the probe names as row names, and the following columns: "chr" (chromosome); "pos" (genomic location of the probe in the genome); and "strand" (this is very important to set up, since the VMRs will only be created based on CpGs on the same strand; if the positions are reported based on a single DNA strand, this should contain either a vector of only "+", "-" or "*" for all of the probes).
#' @param methylation_data A data frame containing M or B values, with samples as columns and probes as rows. Data is expected to have already passed through quality control and cleaning steps.
#' @param cor_threshold Numeric value (0-1) to be used as the median pearson correlation threshold for identifying VMRs (i.e.
#' all VMRs will have a median pairwise probe correlation higher than this threshold).
#' @param var_method A string indicating the metric to use to represent variability in the data set. The options are "mad" (median absolute deviation)
#' or "variance".
#' @param var_distribution A string indicating which probes in the data set should be used to create a variability distribution; the threshold to identify Highly Variable Probes (determined also with the var_threshold_percentile argument) is established based on this distribution. The options 1 is "ultrastable" (a subset of CpGs that are stably methylated/unmethylated across human tissues and developmental states described by [Edgar R., et al.](https://doi.org/10.1186/1756-8935-7-28) in 2014). This option is recommended, especially if you want to compare different populations or tissues, as the threshold value should be comparable. On the other hand, the user can use option 2: "all" (all probes in the data set). The "ultrastable" option is only compatible with Illumina human microarrays. The default is "ultrastable".
#' @param var_threshold_percentile The percentile (0-1) to be used as cutoff to define Highly Variable Probes (which are then grouped into VML). If using the variability of the "ultrastable" probes, we recommend a high threshold (default is 0.99), since these probes are expected to display a very low variation in human tissues. If using the variability of "all" probes, we recommend using a percentile of 0.9 since it captures the top 10% most variable probes, which has been traditionally used in studies. It is important to note that the top 10% most variable probes will capture the same amount of probes in a data set regardless of their overall variability levels, which might differ between tissues or populations.
#' @param max_distance Maximum distance in base pairs allowed for two probes to be grouped into a region. The default is 1000.
#'
#' @return A list with the following elements:
#'  - $var_score_threshold: threshold used to define Highly Variable Probes (mad or variance, depending on the specified choice).
#'  - $highly_variable_probes: a data frame with the probes that passed the variability score threshold imposed by the user, and their variability score (MAD score or variance).
#'  - $VML: a GRanges-like data frame with VMRs (regions composed of two or more contiguous, correlated and proximal Highly Variable Probes), and sVMPs (highly variable probes without neighboring CpGs measured in *max_distance* on the array).
#'
#' @export
#' @examples
#'
#' VML <- RAMEN::findVML(
#'   methylation_data = RAMEN::test_methylation_data,
#'   array_manifest = "IlluminaHumanMethylationEPICv1",
#'   cor_threshold = 0.15,
#'   var_method = "variance",
#'   var_distribution = "ultrastable",
#'   var_threshold_percentile = 0.99,
#'   max_distance = 1000
#' )
#'
findVML <- function(methylation_data,
                    array_manifest,
                    cor_threshold = 0.15,
                    var_method = "variance",
                    var_distribution = "ultrastable",
                    var_threshold_percentile = 0.99,
                    max_distance = 1000) {
  # Check that the array manifest is in the right format
  if (is.data.frame(array_manifest)) {
    if (!all(c("chr", "pos", "strand") %in% colnames(array_manifest))) stop("The array_manifest data frame does not have the required columns. Please provide a manifest with the required columns or provide a string with one of the supported human microarrays ('IlluminaHumanMethylation450k', 'IlluminaHumanMethylationEPICv1','IlluminaHumanMethylationEPICv2')")
    # Check that the array strand is in the format expected by the user
    if (base::length(base::unique(array_manifest$strand)) > 1) warning("The manifest currently has more than one type of strands. Please note that this function is strand sensitive. So, probes in proximal coordinates but different strands on the manifest will not be grouped together. Many array manifests such as the Illumina EPIC one include the PROBE strand, but the position of the actual CpGs (pos) is reported in the same strand; in those cases we recommend setting all of the probes to the same strand.")
    if (var_distribution == "ultrastable") {
      # If the user provides their own manifest and is choosing to use the ultrastable probes, make sure that a good number of them is present in the data set. If not, throw an error
      if (sum(row.names(array_manifest) %in% RAMEN::ultrastable_cpgs) < 100) stop("The var_distribution = 'ultrastable' option is only compatible with Illumina human microarrays at the moment. If you are using a human Illumina microarray please indicate it with their corresponding string, or make sure that it contains a more than 100 ultrastable probes (RAMEN::ultrastable_cpgs). If not, please get the variability threshold based on all the probes in your data set(var_distribution = 'all', var_threshold_percentile = 0.9). ")
    }
  } else if (is.character(array_manifest)) {
    if (!array_manifest %in% c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPICv1", "IlluminaHumanMethylationEPICv2")) stop("The string you provided in array_manifest is not currently supported in RAMEN. Please provide a manifest with the required columns or provide a string with one of the supported human microarrays ('IlluminaHumanMethylation450k', 'IlluminaHumanMethylationEPICv1','IlluminaHumanMethylationEPICv2')")
  } else {
    stop("The array_manifest object is not a data.frame nor a string. Please provide a manifest with the required columns or provide a string with one of the supported human microarrays ('IlluminaHumanMethylation450k', 'IlluminaHumanMethylationEPICv1','IlluminaHumanMethylationEPICv2')")
  }
  #Check that cor_threshold is numeric and between 0 and 1
  if (!is.numeric(cor_threshold) | cor_threshold <= 0 | cor_threshold >= 1) {
    stop("'cor_threshold' must be of type 'numeric' and from 0 to 1")
  }

  # Check that the method choice is correct
  if (var_method == "mad") {
    var_scores <- apply(methylation_data, 1, stats::mad) %>%
      as.data.frame() %>%
      dplyr::rename("var_score" = ".")
  } else if (var_method == "variance") {
    var_scores <- apply(methylation_data, 1, stats::var) %>%
      as.data.frame() %>%
      dplyr::rename("var_score" = ".")
  } else {
    stop("The method must be either 'mad' or 'variance'. Please select one of those options")
  }

  #### Identify highly variable probes ####
  message("Identifying Highly Variable Probes...")
  # Get the variability threshold
  if (var_distribution == "all") {
    var_threshold <- stats::quantile(
      var_scores$var_score,
      var_threshold_percentile
    )
  } else if (var_distribution == "ultrastable") {
    if (is.data.frame(array_manifest)) {
      var_threshold <- stats::quantile(
        var_scores[(row.names(var_scores) %in% RAMEN::ultrastable_cpgs), ], # Subset only ultrastable probes
        var_threshold_percentile
      )
    } else if (array_manifest == "IlluminaHumanMethylationEPICv2") {
      # Get the name of the ultrastable probes in the EPICv2 format
      epicv2_ultrastable_cpgs <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Other |>
        data.frame() |>
        dplyr::filter(Methyl450_Loci %in% RAMEN::ultrastable_cpgs) |>
        tibble::rownames_to_column("epicv2_probes") |>
        dplyr::pull(epicv2_probes)
      var_threshold <- stats::quantile(
        var_scores[(row.names(var_scores) %in% epicv2_ultrastable_cpgs), ],
        var_threshold_percentile
      )
    } else {
      # EPICv1 or 450k (same probe name as the ultrastable probes)
      var_threshold <- stats::quantile(
        var_scores[(row.names(var_scores) %in% RAMEN::ultrastable_cpgs), ], # Subset only ultrastable probes
        var_threshold_percentile
      )
    }
  }
  # Replace the array manifest if the user provided a string with the name of the array
  if (is.character(array_manifest)) {
    if (array_manifest == "IlluminaHumanMethylation450k") {
      manifest <- data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations)
    } else if (array_manifest == "IlluminaHumanMethylationEPICv1") {
      manifest <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations)
    } else if (array_manifest == "IlluminaHumanMethylationEPICv2") {
      manifest <- data.frame(IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations)
    }
  } else {
    manifest <- array_manifest
  }
  # Filter the manifest to remove the probes that have no variability score information because they were not measured/did not pass the QC and are not highly variable
  manifest_hvp <- manifest %>%
    tibble::rownames_to_column(var = "TargetID") %>%
    dplyr::select(c(TargetID, chr, pos, strand)) %>%
    dplyr::filter(
      !is.na(pos), # Remove probes with no map info
      TargetID %in% row.names(var_scores %>%
                                dplyr::filter(var_score >= var_threshold))
      ) %>% # Remove probes that have no methylation information in the processed data and are not highly variable
    dplyr::left_join(
      var_scores %>% # Add variability information
        tibble::rownames_to_column(var = "TargetID"),
      by = "TargetID"
    ) %>%
    dplyr::arrange(chr) %>% # important step for using Rle later when constructing the GenomicRanges object!
    as.data.frame()
  rownames(manifest_hvp) <- manifest_hvp$TargetID
  if (is.factor(manifest_hvp$chr)) manifest_hvp <- manifest_hvp %>% dplyr::mutate(chr = droplevels(chr))

  #### Identify sparse Variable Methylated Probes####
  message("Identifying sparse Variable Methylated Probes")
  full_manifest <- manifest %>%
    tibble::rownames_to_column(var = "TargetID") %>%
    dplyr::select(c(TargetID, chr, pos, strand)) %>%
    dplyr::filter(
      !is.na(pos), # Remove probes with no map info
      TargetID %in% row.names(var_scores)
    ) %>% # keep only the probes where we have methylation information
    dplyr::arrange(chr) %>% # important step for using Rle later when constructing the GenomicRanges object!
    as.data.frame()
  rownames(full_manifest) <- full_manifest$TargetID
  if (is.factor(full_manifest$chr)) full_manifest <- full_manifest %>% dplyr::mutate(chr = droplevels(chr))

  # Convert the full manifest to a GenomicRanges object
  seqnames_full_manifest_gr <- table(full_manifest$chr)
  full_manifest_gr <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(names(seqnames_full_manifest_gr), as.numeric(seqnames_full_manifest_gr)), # Number of chromosome; as.numeric to convert from table to numeric vector
    ranges = IRanges::IRanges(full_manifest$pos,
      end = full_manifest$pos,
      names = full_manifest$TargetID
    ),
    strand = S4Vectors::Rle(
      rle(as.character(full_manifest$strand))$values,
      rle(as.character(full_manifest$strand))$lengths
    )
  )

  #### Group the probes into regions to detect sVMPs####
  regions_full_manifest <- GenomicRanges::reduce(full_manifest_gr, with.revmap = TRUE, min.gapwidth = max_distance)
  # Add the number of probes in each region
  S4Vectors::mcols(regions_full_manifest)$n_probes <- sapply(S4Vectors::mcols(regions_full_manifest)$revmap, length)
  # Substitute revmap with the name of the probes in each region
  S4Vectors::mcols(regions_full_manifest)$probes <- sapply(S4Vectors::mcols(regions_full_manifest)$revmap, map_revmap_names, full_manifest)
  # Remove revmap mcol
  S4Vectors::mcols(regions_full_manifest)$revmap <- NULL
  # Keep elements with only one probe
  lonely_probes <- regions_full_manifest[(GenomicRanges::elementMetadata(regions_full_manifest)[, "n_probes"] <= 1)] %>%
    as.data.frame() %>%
    dplyr::pull(probes) %>%
    unlist()

  #### Identify VMRs####
  message("Identifying Variable Methylated Regions...")
  # convert the highly variable probes data frame to a GenomicRanges object
  seqnames_gr <- table(manifest_hvp$chr)
  gr <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(names(seqnames_gr), as.numeric(seqnames_gr)), # Number of chromosome; as.numeric to convert from table to numeric vector
    ranges = IRanges::IRanges(manifest_hvp$pos,
      end = manifest_hvp$pos,
      names = manifest_hvp$TargetID
    ),
    strand = S4Vectors::Rle(
      rle(as.character(manifest_hvp$strand))$values,
      rle(as.character(manifest_hvp$strand))$lengths
    ),
    var_score = manifest_hvp$var_score
  ) # Metadata

  # Create the regions
  candidate_VMRs <- GenomicRanges::reduce(gr, with.revmap = TRUE, min.gapwidth = max_distance)
  # Add the number of probes in each region
  S4Vectors::mcols(candidate_VMRs)$n_VMPs <- sapply(S4Vectors::mcols(candidate_VMRs)$revmap, length)
  # Substitute revmap with the name of the probes in each VMR
  S4Vectors::mcols(candidate_VMRs)$probes <- sapply(S4Vectors::mcols(candidate_VMRs)$revmap, map_revmap_names, manifest_hvp)
  # Remove revmap mcol
  S4Vectors::mcols(candidate_VMRs)$revmap <- NULL

  ### Capture canonical VMRs ###
  message("Applying correlation filter to Variable Methylated Regions...")
  VMRs <- candidate_VMRs[(GenomicRanges::elementMetadata(candidate_VMRs)[, "n_VMPs"] > 1)] %>%
    data.frame() # Convert the GR to a data frame so that I can use medCorVMR()
  ### Check for correlation between probes only if we have VMRs
  if (nrow(VMRs) > 0) {
    VMRs <- VMRs %>%
      medCorVMR(VMR_df = ., methylation_data = methylation_data) %>% # Compute the median correlation of each region
      dplyr::filter(median_correlation > cor_threshold) %>% # Remove VMRs whose CpGs are not correlated
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) # Create a GR object again
  } else {
    warning("No canonical VMRs were found in this data set")
  }

  ### Capture non-canonical VMRs ###
  sVMPs <- candidate_VMRs[(GenomicRanges::elementMetadata(candidate_VMRs)[, "probes"] %in% lonely_probes)] # Select the lonely probes
  GenomicRanges::mcols(sVMPs)$median_correlation <- rep(NA, nrow(GenomicRanges::mcols(sVMPs))) # Add a column of NAs under the name of median_correlation to match the strict_VMRs

  return(list(
    var_score_threshold = var_threshold,
    highly_variable_probes = var_scores %>%
      tibble::rownames_to_column(var = "TargetID") %>%
      dplyr::filter(TargetID %in% manifest_hvp$TargetID),
    VML = data.frame(VMRs) %>%
      rbind(data.frame(sVMPs)) %>%
      dplyr::mutate(
        type = ifelse(n_VMPs > 1, "VMR", "sVMP"),
        VML_index = paste("VML", as.character(dplyr::row_number()), sep = "")
      ) %>%
      dplyr::select(VML_index, type, seqnames, start, end, width, strand, probes, n_VMPs, median_correlation)
  ))
}
