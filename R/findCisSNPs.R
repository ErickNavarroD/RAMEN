#' Find cis SNPs around a set of Variable Methylated Loci (VML)
#'
#' Identification of genotyped Single Nucleotide Polymorphisms (SNPs) close to
#' each VML using a distance threshold.
#'
#' **Important**: please make sure that the positions of the VML data frame and
#' the ones in the genotype information are from the same genome build.
#'
#' @param VML_df A GRanges-like data frame (i.e. the same columns as a GRanges
#' object converted to a data frame). Must contain the following columns:
#' "seqnames", "start", "end". These columns are present automatically when
#' doing the object conversion and correspond to the chromosome number, and
#' range of the region.
#' @param genotype_information A data frame with information about genotyped
#' sites of interest. It must contain the following columns: "CHROM"
#' (chromosome number), "POS" (Genomic basepair position of the SNP (must be an
#' integer), and "ID" (SNP ID). The nomenclature of CHROM must match with the
#' one used in the VML_df seqnames column (i.e., if VML_df$seqnames uses 1, 2,
#' 3, X, Y or Chr1, Chr2, Chr3, ChrX, ChrY, etc. as chromosome number, the
#' genotype_information$CHROM values must be encoded in the same way).
#' @param distance The distance threshold in basepairs to be used to identify
#' cis SNPs. Default is 1 Mb.
#'
#' @return The same VML data frame (a data frame compatible with GRanges
#' conversion) with the following new columns:
#'  - The cis SNPs identified for each VML and the number of SNPs surrounding
#'  each VML in the specified window
#' @export
#' @examples
#' ## Find VML in test data
#' VML <- RAMEN::findVML(
#'    methylation_data = RAMEN::test_methylation_data,
#'    array_manifest = "IlluminaHumanMethylationEPICv1",
#'    cor_threshold = 0,
#'    var_method = "variance",
#'    var_distribution = "ultrastable",
#'    var_threshold_percentile = 0.99,
#'    max_distance = 1000
#'    )
#' ## Find cis SNPs around VML
#' VML_with_cis_snps <- RAMEN::findCisSNPs(
#'   VML_df = VML$VML,
#'   genotype_information = RAMEN::test_genotype_information,
#'   distance = 1e6
#'   )
#'
findCisSNPs <- function(VML_df, genotype_information, distance = 1e6) {
  CHROM <- NULL
  # Check arguments
  if (!is.data.frame(VML_df)) stop("Please make sure the VML_df object is a data frame.")
  if (!is.data.frame(genotype_information)) stop("Please make sure the genotype_information object is a data frame.")
  if (!all(c("seqnames", "start", "end") %in% colnames(VML_df))) stop("Please make sure the VML_df object has the required columns with the appropiate names (check documentation for further information)")
  if (!all(c("CHROM", "POS", "ID") %in% colnames(genotype_information))) stop("Please make sure the genotype_information object has the required columns with the appropiate names (check documentation for further information)")
  message("Reminder: please make sure that the positions of the VML data frame and the ones in the genotype information are from the same genome build.")
  # Convert VML and snp data into a GenomicRanges object
  VML_gr <- GenomicRanges::makeGRangesFromDataFrame(VML_df, keep.extra.columns = TRUE)
  genotype_information <- genotype_information %>%
    dplyr::arrange(CHROM) # important step for using Rle later when constructing the GenomicRanges object!
  seqnames_gr <- table(genotype_information$CHROM)
  genot_gr <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(names(seqnames_gr), as.numeric(seqnames_gr)), # Number of chromosome; as.numeric to convert from table to numeric vector
    ranges = IRanges::IRanges(genotype_information$POS,
      end = genotype_information$POS,
      names = genotype_information$ID
    )
  )
  # Extend each VML 1 Mb up and downstream
  VML_extended <- VML_gr + distance

  VML_df_with_cisSNPs <- VML_df
  if (!"VML_index" %in% colnames(VML_df_with_cisSNPs)) { # Add a VML index to each region if not already existing
    VML_df_with_cisSNPs <- VML_df_with_cisSNPs %>%
      dplyr::mutate(VML_index = paste("VML", as.character(dplyr::row_number()), sep = ""))
  }

  #### Get the number of overlaps per extended VML ####
  VML_df_with_cisSNPs$surrounding_SNPs <- GenomicRanges::countOverlaps(VML_extended, genot_gr)

  #### Identify the SNPs that are present in each VML ####
  snps_per_vml_find <- GenomicRanges::findOverlaps(VML_extended, genot_gr, select = "all")
  rownames(genotype_information) <- genotype_information$ID
  VML_df_with_cisSNPs <- VML_df_with_cisSNPs %>%
    dplyr::mutate(SNP = lapply(snps_per_vml_find, map_revmap_names, genotype_information))

  return(VML_df_with_cisSNPs)
}
