#' Summarize the methylation states of Variable Methylated Loci (VML)
#'
#' This function computes a representative methylation score for each Variable
#' Methylated Locus (VML) in a dataset. It returns a data frame with the median
#' methylation of each region per individual. For each VML in a dataset, returns
#' a  with the median methylation of that region (columns) per individual
#' (rows) as representative score.
#'
#' This function supports parallel computing for increased speed. To do so, you
#' have to set the parallel backend in your R session BEFORE running the
#' function (e.g., *doParallel::registerDoParallel(4)*). After that, the
#' function can be run as usual.
#'
#' @inheritParams medCorVMR
#'
#' @return A data frame with samples as rows, and VML as columns. The value
#' inside each cell corresponds to the summarized methylation value of said VML
#' in the corresponding individual. The column names correspond to the VML_index.
#'
#' @importFrom foreach %dopar%
#' @export
#'
#' @examples
#' ## Find VML in test data
#' # Set the parallel backend to use 2 workers
#' doParallel::registerDoParallel(2)
#' VML <- findVML(
#'   methylation_data = test_methylation_data,
#'   array_manifest = "IlluminaHumanMethylationEPICv1",
#'   cor_threshold = 0,
#'   var_method = "variance",
#'   var_distribution = "ultrastable",
#'   var_threshold_percentile = 0.99,
#'   max_distance = 1000
#' )
#'
#' ## Summarize methylation states of the found VML
#' summarized_VML <- summarizeVML(
#'   # Use only 5 for demonstration purposes
#'   VML = VML$VML[1:5, ],
#'   methylation_data = test_methylation_data
#' )
#'
summarizeVML <- function(VML,
                         methylation_data) {
  #### Input checks ####
  argument_check(VML, "GRanges")
  argument_check(methylation_data, "data.frame")
  # Add a VML index to each region if not already existing
  if (!"VML_index" %in% colnames(S4Vectors::mcols(VML))) {
    S4Vectors::mcols(VML)$VML_index = paste("VML", 1:length(VML), sep = "")
  }

  if (!all(unique(unlist(VML$probes)) %in% rownames(methylation_data))) {
    warning(paste("Some probes listed in the VML data frame are not found in",
    "the methylation data. Please check that all probes listed in the 'probes'",
    "column of the VML data frame are present in the row names of the",
    "methylation data frame to avoid having NAs."))
  }
  # Check that probes is a list.
  if (!is.list(VML$probes)) {
    stop("Please make sure the 'probes' column in the VML object is a list")
  }

  VML_index <- i <- NULL # To avoid R CMD check notes

  #### Summarize VML ####
  summarized_VML <- foreach::foreach(i = VML$VML_index,
                                     .combine = "cbind") %dopar% {
    probes <- VML[VML$VML_index == i]$probes |>
      unlist()
    subset_meth <-  methylation_data[probes, , drop = FALSE]
    median <- apply(subset_meth, 2, median, na.rm = TRUE)
    matrix(median, ncol = 1, dimnames = list(NULL, i))
                                     }
  # Add ID names
  rownames(summarized_VML) = colnames(methylation_data)
  return(data.frame(summarized_VML))
}
