# Compute the median probe methylation pearson correlation for each Variable Methylated Region (VMR).

This function will take a GRanges object converted into a data frame,
where each row corresponds to a Variable Methylated Region. Then, it
computes the pairwise correlation of the probes of each VMR and reports
its median pairwise probe correlation.

## Usage

``` r
medCorVMR(VML, methylation_data)
```

## Arguments

- VML:

  GRanges object. Must contain a metadata column named "probes", where
  each element contains a vector with the probes constituting the VML.

- methylation_data:

  A data frame containing M or B values, with samples as columns and
  probes as rows. Data is expected to have already passed through
  quality control and cleaning steps. Rows must be the CpG probe IDs.

## Value

A GRanges object like VML with an extra column per region containing the
median pairwise correlation.

## Details

This function supports parallel computing for increased speed. To do so,
you have to set the parallel backend in your R session before running
the function (e.g., *doParallel::registerDoParallel(4)*)). After that,
the function can be run as usual. It is recommended to also set
options(future.globals.maxSize= +Inf).

## Examples

``` r
# Evaluate sequentially
foreach::registerDoSEQ()
# Create a VML object
VML <- GenomicRanges::GRanges(seqnames = c("chr21", "chr21"),
      ranges = IRanges::IRanges(start = c(10861376, 10862171),
                       end = c(10862507, 10883548)),
      probes = I(list(
        c("cg15043638", "cg18287590", "cg17975851"),
        c("cg13893907", "cg17035109", "cg06187584")))
      )

# Compute median correlation for each VMR
medCorVMR(VML = VML, methylation_data = RAMEN::test_methylation_data)
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames            ranges strand |                           probes
#>          <Rle>         <IRanges>  <Rle> |                           <list>
#>   [1]    chr21 10861376-10862507      * | cg15043638,cg18287590,cg17975851
#>   [2]    chr21 10862171-10883548      * | cg13893907,cg17035109,cg06187584
#>       median_correlation
#>                <numeric>
#>   [1]           0.727516
#>   [2]           0.746814
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
