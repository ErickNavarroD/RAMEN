#' Methylation data matrix example
#'
#' Simulated M values of the 3000 probes selected in *test_array_manifest* for 30 individuals. Values were converted from beta values, which were drawn from a bimodal Beta distribution.
#'
#' @format ## `test_methylation_data`
#' A data frame with 3,000 rows and 30 columns:
#' \describe{
#'   \item{*rownames*}{Probe IDs (column TargetID in the EPIC array)}
#'   \item{ID1:30}{DNAme profile of individuals 1 to 30; column names correspond to individual IDs}
#'   ...
#' }
"test_methylation_data"
