#' Find cis SNPs around a set of VMRs
#'
#' Identification of genotyped SNPs close to each VMR using a distance threshold. Important: please make sure that the
#' positions of the VMR data frame and the ones in the genotype information are from the same genome build.
#'
#' @param VMRs_df A data frame converted from a Genomic Ranges object. Must contain the following columns:
#' "seqnames", "start", "end".
#' @param genotype_information A data frame with information about genotyped sites of interest. It must contain the following
#' columns: "Chr" - Number of chromosome, "Position" - Genomic position of the chromosome (must contain values of class int),
#' and "Name" - ID of the SNP.
#' @param distance The distance threshold to be used to identify cis SNPs
#'
#' @return a VMR_df: a VMR_df object but now with a new column including the cis SNPs identified for each VMR, and the number
#' of SNPs surrounding each VMR in the specified window
#' @export

findCisSNPs = function(VMRs_df, genotype_information, distance = 1e6){
  #Convert VMR and snp data into a GenomicRanges object
  VMRs_gr = GenomicRanges::makeGRangesFromDataFrame(VMRs_df, keep.extra.columns = TRUE)
  genotype_information = genotype_information %>% dplyr::arrange(Chr) #important step for using Rle later when constructing the GenomicRanges object!
  seqnames_gr = table(genotype_information$Chr)
  genot_gr = GenomicRanges::GRanges(
    seqnames =  S4Vectors::Rle(names(seqnames_gr), as.numeric(seqnames_gr)), #Number of chromosome; as.numeric to convert from table to numeric vector
    ranges = IRanges::IRanges(genotype_information$Position, end = genotype_information$Position ,
                              names = genotype_information$Name))
  #Extend each VMR 1 Mb up and downstream
  VMRs_extended = VMRs_gr + distance

  VMRs_df_with_cisSNPs = VMRs_df
  if(!"VMR_index" %in% colnames(VMRs_df_with_cisSNPs)){ # Add a VMR index to each region if not already existing
    VMRs_df_with_cisSNPs = VMRs_df_with_cisSNPs %>%
      tibble::rownames_to_column(var = "VMR_index")
  }

  #### Get the number of overlaps per extended VMR ####
  VMRs_df_with_cisSNPs$surrounding_SNPs =  GenomicRanges::countOverlaps(VMRs_extended, genot_gr)

  ####Identify the SNPs that are present in each VMR ####
  snps_per_vmr_find =  GenomicRanges::findOverlaps(VMRs_extended, genot_gr, select = "all")
  rownames(genotype_information) = genotype_information$Name
  VMRs_df_with_cisSNPs = VMRs_df_with_cisSNPs %>%
    dplyr::mutate(SNP = sapply(snps_per_vmr_find, map_revmap_names, genotype_information))

  return(VMRs_df_with_cisSNPs)
}
