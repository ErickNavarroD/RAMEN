#' Find cis SNPs around a set of VMRs
#'
#' Identification of genotyped SNPs close to each VMR using a distance threshold. Important: please make sure that the
#' positions of the VMR data frame and the ones in the genotype information are from the same genome build.
#'
#' @param VMRs_df A data frame converted from a Genomic Ranges object. Must contain the following columns:
#' "probes", containing a list where each element contains a vector with the probes constituting a  VMR
#' @param genotype_information A data frame with information about genotyped sites of interest. It must contain the following
#' columns: "Chr" - Number of chromosome, "Position" - Genomic position of the chromosome, and "Name" - ID of the SNP.
#' @param distance The distance threshold to be used to identify cis SNPs
#'
#' @return a list with the following elements:
#' - A GenomicRanges::countOverlaps() output, which is an integer vector containing the tabulated query overlap hits.
#'  In other words, it is a data frame with two columns; the first one is the VMR number ID (which corresponds
#'  to the order of the VMRs in the VMR_df object), and a second one with the number of cis SNPs that were found. This
#'  is a useful object for ploting cisSNP distributions.
#' - VMR_df: a VMR_df object but now with a new column including the cis SNPs identified for each VMR.
#' @export

findCisSNPs = function(VMRs_df, genotype_information, distance = 1e6){
  #Convert VMR and snp data into a GenomicRanges object
  VMRs_gr = GenomicRanges::makeGRangesFromDataFrame(VMRs_df, keep.extra.columns = TRUE)
  seqnames_gr = table(genotype_information$Chr)
  genot_gr = GenomicRanges::GRanges(
    seqnames =  S4Vectors::Rle(names(seqnames_gr), as.numeric(seqnames_gr)), #Number of chromosome
    ranges = IRanges::IRanges(genotype_information$Position, end = genotype_information$Position ,
                              names = genotype_information$Name))
  #Extend each VMR 1 Mb up and downstream
  VMRs_extended = VMRs_gr + distance
  #### Get the number of overlaps per extended VMR ####
  snps_per_vmr_counts = data.frame(VMR_index = c(1:length(VMRs_extended)),
                                   number_of_SNPs =  GenomicRanges::countOverlaps(VMRs_extended, genot_gr))

  ####Identify the SNPs that are present in each VMR ####
  ### this is just for the purpose of the preliminary analysis
  #snps_per_vmr_find = findOverlaps(VMRs_extended, genot_gr, select = "all")
  snps_per_vmr_find =  GenomicRanges::findOverlaps(VMRs_extended, genot_gr, select = "first")
  rownames(genotype_information) = genotype_information$Name
  VMRs_df_with_cisSNPs = VMRs_df %>%
    tibble::rownames_to_column(var = "VMR_index") %>%
    dplyr::mutate(SNP = sapply(snps_per_vmr_find, map_revmap_names, genotype_information))

  return(list(snps_per_vmr_counts = snps_per_vmr_counts,
              VMRs_df_with_cisSNPs = VMRs_df_with_cisSNPs))
}
