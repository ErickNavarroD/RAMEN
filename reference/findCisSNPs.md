# Find cis SNPs around a set of Variable Methylated Regions (VMRs)

Identification of genotyped Single Nucleotide Polymorphisms (SNPs) close
to each VMR using a distance threshold.

## Usage

``` r
findCisSNPs(VMRs_df, genotype_information, distance = 1e+06)
```

## Arguments

- VMRs_df:

  A GRanges object converted to a data frame. Must contain the following
  columns: "seqnames", "start", "end". These columns are present
  automatically when doing the object conversion and correspond to the
  chromosome number, and range of the region.

- genotype_information:

  A data frame with information about genotyped sites of interest. It
  must contain the following columns: "CHROM" - chromosome number,
  "POS" - Genomic basepair position of SNP in the corresponding
  chromosome (must contain values of class int), and "ID" - SNP ID. The
  nomenclature of CHROM must match with the one used in the VMRs_df
  seqnames column (i.e., if VMRs_df\$seqnames uses 1, 2, 3, X, Y or
  Chr1, Chr2, Chr3, ChrX, ChrY, etc. as chromosome number, the
  genotype_information\$CHROM values must be encoded in the same way).

- distance:

  The distance threshold to be used to identify cis SNPs. Default is 1
  Mb.

## Value

A VMR_df object (a data frame compatible with GRanges conversion) with
the following new columns:

- The cis SNPs identified for each VMR, the number of SNPs surrounding
  each VMR in the specified window

- VMR_index, which is created if not already existing based on the
  rownames of the VMR_df.

## Details

**Important**: please make sure that the positions of the VMR data frame
and the ones in the genotype information are from the same genome build.
