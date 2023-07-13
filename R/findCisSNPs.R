#' Find cis SNPs around a set of Variable Methylated Regions (VMRs)
#'
#' Identification of genotyped Single Nucleotide Polymorphisms (SNPs) close to each VMR using a distance threshold.
#'
#' **Important**: please make sure that the positions of the VMR data frame and the ones in the genotype information are from the same genome build.
#'
#' @param VMRs_df A GRanges object converted to a data frame. Must contain the following columns:
#' "seqnames", "start", "end". These columns are present automatically when doing the object conversion and correspond to the chromosome number, and range of the region.
#' @param genotype_information A data frame with information about genotyped sites of interest. It must contain the following
#' columns: "CHROM" - chromosome number, "POS" - Genomic basepair position of SNP in the corresponding
#' chromosome (must contain values of class int), and "ID" - SNP ID. The nomenclature of CHROM must match with the one used in the VMRs_df seqnames column (i.e., if VMRs_df$seqnames uses 1, 2, 3, X, Y or Chr1, Chr2, Chr3, ChrX, ChrY, etc. as chromosome number, the genotype_information$CHROM values must be encoded in the same way).
#' @param distance The distance threshold to be used to identify cis SNPs. Default is 1 Mb.
#'
#' @return A VMR_df object (a data frame compatible with GRanges conversion) with the following new columns:
#'  - The cis SNPs identified for each VMR, the number of SNPs surrounding each VMR in the specified window
#'  - VMR_index, which is created if not already existing based on the rownames of the VMR_df.
#' @export

findCisSNPs = function(VMRs_df, genotype_information, distance = 1e6){
  #Check arguments
  if(!all(c("seqnames","start","end") %in% colnames(VMRs_df))) stop("Please make sure the VMRs_df object has the required columns with the appropiate names (check documentation for further information)")
  if(!all(c("CHROM","POS","ID") %in% colnames(genotype_information))) stop("Please make sure the genotype_information object has the required columns with the appropiate names (check documentation for further information)")
  message("Important: please make sure that the positions of the VMR data frame and the ones in the genotype information are from the same genome build.")
  #Convert VMR and snp data into a GenomicRanges object
  VMRs_gr = GenomicRanges::makeGRangesFromDataFrame(VMRs_df, keep.extra.columns = TRUE)
  genotype_information = genotype_information %>% dplyr::arrange(CHROM) #important step for using Rle later when constructing the GenomicRanges object!
  seqnames_gr = table(genotype_information$CHROM)
  genot_gr = GenomicRanges::GRanges(
    seqnames =  S4Vectors::Rle(names(seqnames_gr), as.numeric(seqnames_gr)), #Number of chromosome; as.numeric to convert from table to numeric vector
    ranges = IRanges::IRanges(genotype_information$POS, end = genotype_information$POS ,
                              names = genotype_information$ID))
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
  rownames(genotype_information) = genotype_information$ID
  VMRs_df_with_cisSNPs = VMRs_df_with_cisSNPs %>%
    dplyr::mutate(SNP = sapply(snps_per_vmr_find, map_revmap_names, genotype_information))

  return(VMRs_df_with_cisSNPs)
}
