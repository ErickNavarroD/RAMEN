# Find cis SNPs around a set of Variable Methylated Loci (VML)

Identification of genotyped Single Nucleotide Polymorphisms (SNPs) close
to each VML using a distance threshold.

## Usage

``` r
findCisSNPs(VML, genotype_information, distance = 1e+06)
```

## Arguments

- VML:

  GRanges object. Must contain a metadata column named "probes", where
  each element contains a vector with the probes constituting the VML.

- genotype_information:

  A data frame with information about genotyped sites of interest. It
  must contain the following columns: "CHROM" (chromosome number), "POS"
  (Genomic basepair position of the SNP (must be an integer), and "ID"
  (SNP ID). The nomenclature of CHROM must match with the one used in
  the VML seqnames column (i.e., if VML uses 1, 2, 3, X, Y or Chr1,
  Chr2, Chr3, ChrX, ChrY, etc. as chromosome number, the
  genotype_information\$CHROM values must be encoded in the same way).

- distance:

  The distance threshold in basepairs to be used to identify cis SNPs.
  Default is 1 Mb.

## Value

The same VML object with new metadata columns indicating the cis SNPs
identified for each VML and the number of SNPs surrounding each VML in
the specified window

## Details

For DNAme data, previous studies (e.g. Gibbs, et al., 2010, McClay, et
al., 2015) have found that SNPs are more likely to associate with DNAme
levels the closer they are to the CpG. We recommend to include SNPs
within 500 kb to 1 Mb to capture SNPs with a high potential to associate
with DNAme.

**Important**: please make sure that the positions of the VML data frame
and the ones in the genotype information are from the same genome build.

## Examples

``` r
## Find VML in test data
VML <- findVML(
  methylation_data = RAMEN::test_methylation_data,
  array_manifest = "IlluminaHumanMethylationEPICv1",
  cor_threshold = 0,
  var_method = "variance",
  var_distribution = "ultrastable",
  var_threshold_percentile = 0.99,
  max_distance = 1000
)
#> Identifying Highly Variable Probes...
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
#> Identifying sparse Variable Methylated Probes
#> Identifying Variable Methylated Regions...
#> Applying correlation filter to Variable Methylated Regions...
#> Warning: executing %dopar% sequentially: no parallel backend registered
## Find cis SNPs around VML
# Use only 5 for demonstration purposes
VML_with_cis_snps <- findCisSNPs(
  VML = VML$VML[1:5, ],
  genotype_information = RAMEN::test_genotype_information,
  distance = 1e6
)
#> Reminder: please make sure that the positions of the VML data frame and the ones in the genotype information are from the same genome build.
```
