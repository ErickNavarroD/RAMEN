#' Array manifest example data set
#'
#' A subset of data from Illumina's EPIC array manifest (first 3,000 probes of the chromosome 21).
#'
#' @format ## `test_array_manifest`
#' A data frame with 3,000 rows and 3 columns:
#' \describe{
#'   \item{*rownames*}{Probe ID - for storage reasons, this variable was stored as row names, but rownames have to be converted to a new column called "TargetID" prior to its use in RAMEN.}
#'   \item{MAPINFO}{Probe genomic position (h19)}
#'   \item{CHR}{Chromosome}
#'   \item{STRAND}{Strand}
#'   ...
#' }
#' @source <https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip>
"test_array_manifest"
