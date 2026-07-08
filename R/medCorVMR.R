#' Compute the median probe methylation pearson correlation for each Variable
#' Methylated Region (VMR).
#'
#' This function will take a GRanges object converted into a data frame, where
#' each row corresponds to a Variable Methylated Region. Then, it computes the
#' pairwise correlation of the probes of each VMR and reports its median
#' pairwise probe correlation.
#'
#' This function supports parallel computing for increased speed. To do so, you
#' have to set the parallel backend in your R session before running the
#' function (e.g., *doParallel::registerDoParallel(4)*)). After that, the
#' function can be run as usual. It is recommended to also set
#' options(future.globals.maxSize= +Inf).
#'
#' @param VML GRanges object. Must contain a metadata column named "probes",
#' where each element contains a vector with the probes constituting the
#' VML.
#' @inheritParams findVML
#' @return A GRanges object like VML with an extra column per region containing
#' the median pairwise correlation.
#'
#' @importFrom foreach %dopar%
#' @export
#'
#' @examples
#' # Set the parallel backend to use 2 workers
#' doParallel::registerDoParallel(2)
#' # Create a VML object
#' VML <- GenomicRanges::GRanges(seqnames = c("chr21", "chr21"),
#'       ranges = IRanges::IRanges(start = c(10861376, 10862171),
#'                        end = c(10862507, 10883548)),
#'       probes = I(list(
#'         c("cg15043638", "cg18287590", "cg17975851"),
#'         c("cg13893907", "cg17035109", "cg06187584")))
#'       )
#'
#' # Compute median correlation for each VMR
#' medCorVMR(VML = VML, methylation_data = RAMEN::test_methylation_data)
#'
medCorVMR <- function(VML, methylation_data) {
  argument_check(VML, "GRanges")
  if (!"probes" %in% colnames(mcols(VML))) {
    stop("Please make sure the VML object has the 'probes' column.")
  }
  argument_check(methylation_data, "data.frame")
  # generate a list where each element will contain a vector with the probes
  # present in one VMR
  VMR_probes <- S4Vectors::mcols(VML)$probes
  # Compute correlations
  i <- NULL # Bind variable to the environment
  median_correlation <- foreach::foreach(
    i = seq_along(VMR_probes), # For each VMR
    .combine = "c" # Combine outputs in a vector
  ) %dopar% {
    # If the VMR has one probe
    if (length(VMR_probes[[i]]) == 1) {
      NA
    } else {
      VMR_correlation <- c()
      # For each probe except the last one
      for (probe_x_i in 1:(length(VMR_probes[[i]]) - 1)) {
        primary_probe <- VMR_probes[[i]][probe_x_i]
        # compute the pairwise correlation with the downstream probes
        for (probe_y_i in (probe_x_i + 1):length(VMR_probes[[i]])) {
          secondary_probe <- VMR_probes[[i]][probe_y_i]
          VMR_correlation <- c(
            VMR_correlation,
            # unlist added to make the subset df a vector
            stats::cor(unlist(methylation_data[primary_probe, ]),
              unlist(methylation_data[secondary_probe, ]),
              method = "pearson"
            )
          )
        }
      }
      median_correlation <- stats::median(VMR_correlation)
    }
  }

  S4Vectors::mcols(VML)$median_correlation <- median_correlation
  return(VML)
}
