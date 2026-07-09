#' Find cis SNPs around a set of Variable Methylated Loci (VML)
#'
#' Identification of genotyped Single Nucleotide Polymorphisms (SNPs) close to
#' each VML using a distance threshold.
#'
#' **Important**: please make sure that the positions of the VML data frame and
#' the ones in the genotype information are from the same genome build.
#'
#' @inheritParams medCorVMR
#' @param genotype_information A data frame with information about genotyped
#' sites of interest. It must contain the following columns: "CHROM"
#' (chromosome number), "POS" (Genomic basepair position of the SNP (must be an
#' integer), and "ID" (SNP ID). The nomenclature of CHROM must match with the
#' one used in the VML seqnames column (i.e., if VML uses 1, 2,
#' 3, X, Y or Chr1, Chr2, Chr3, ChrX, ChrY, etc. as chromosome number, the
#' genotype_information$CHROM values must be encoded in the same way).
#' @param distance The distance threshold in basepairs to be used to identify
#' cis SNPs. Default is 1 Mb.
#'
#' @return The same VML object with new metadata columns indicating the cis SNPs
#'  identified for each VML and the number of SNPs surrounding each VML in the
#'  specified window
#' @export
#' @examples
#' ## Find VML in test data
#' VML <- RAMEN::findVML(
#'   methylation_data = RAMEN::test_methylation_data,
#'   array_manifest = "IlluminaHumanMethylationEPICv1",
#'   cor_threshold = 0,
#'   var_method = "variance",
#'   var_distribution = "ultrastable",
#'   var_threshold_percentile = 0.99,
#'   max_distance = 1000
#' )
#' ## Find cis SNPs around VML
#' # Use only 5 for demonstration purposes
#' VML_with_cis_snps <- RAMEN::findCisSNPs(
#'   VML = VML$VML[1:5, ],
#'   genotype_information = RAMEN::test_genotype_information,
#'   distance = 1e6
#' )
#'
findCisSNPs <- function(VML, genotype_information, distance = 1e6) {
  CHROM <- NULL
  #### Check arguments ####
  argument_check(VML, "GRanges")
  if (!"probes" %in% colnames(mcols(VML))) {
    stop("Please make sure the VML object has the 'probes' column.")
  }
  argument_check(genotype_information, "data.frame")
  columns_exist(genotype_information, c("CHROM", "POS", "ID"))
  message(paste("Reminder: please make sure that the positions of the VML data",
  "frame and the ones in the genotype information are from the same genome",
  "build."))
  argument_check(distance, "numeric")

  #### Extend each VML 1 Mb up and downstream ####
  #Convert genotype information into GRanges object
  genotype_information <- genotype_information |>
    # important step for using Rle later when constructing the GenomicRanges
    # object
    dplyr::arrange(CHROM)
  seqnames_gr <- table(genotype_information$CHROM)
  genot_gr <- GenomicRanges::GRanges(
    # Number of chromosome; as.numeric to convert from table to numeric vector
    seqnames = S4Vectors::Rle(names(seqnames_gr), as.numeric(seqnames_gr)),
    ranges = IRanges::IRanges(genotype_information$POS,
      end = genotype_information$POS,
      names = genotype_information$ID
    )
  )
  # Extend window
  VML_extended <- VML + distance

  # Add a VML index to each region if not already existing
  if (!"VML_index" %in% colnames(S4Vectors::mcols(VML))) {
    S4Vectors::mcols(VML)$VML_index = paste("VML", 1:length(VML), sep = "")
  }

  #### Identify the SNPs that are present in each VML ####
  # Add the number of surrounding SNPs
  VML$surrounding_SNPs <- GenomicRanges::countOverlaps(VML_extended,
                                                       genot_gr)
  snps_per_vml_find <- GenomicRanges::findOverlaps(VML_extended,
                                                   genot_gr,
                                                   select = "all")
  rownames(genotype_information) <- genotype_information$ID
  # Add the IDs of the surrounding SNPs
  VML$SNP <- lapply(snps_per_vml_find,
                    map_revmap_names,
                    genotype_information)

  return(VML)
}
