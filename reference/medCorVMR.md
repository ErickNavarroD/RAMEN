# Compute the median probe methylation pearson correlation for each Variable Methylated Region (VMR).

This function will take a GRanges object converted into a data frame,
where each row corresponds to a Variable Methylated Region. Then, it
computes the pairwise correlation of the probes of each VMR and reports
its median pairwise probe correlation.

## Usage

``` r
medCorVMR(VMR_df, methylation_data)
```

## Arguments

- VMR_df:

  GRanges object converted to a data frame. Must contain the following
  columns: "seqnames", "start", "end" (all of which are produced
  automatically when doing the object conversion) and "probes"
  (containing a list in which each element contains a vector with the
  probes constituting the VMR).

- methylation_data:

  A data frame containing M or B values, with samples as columns and
  probes as rows. Data is expected to have already passed through
  quality control and cleaning steps.

## Value

A data frame like VMR_df with an extra column per region containing the
median pairwise correlation.

## Details

This function supports parallel computing for increased speed. To do so,
you have to set the parallel backend in your R session before running
the function (e.g., *doParallel::registerDoParallel(4)*)). After that,
the function can be run as usual. It is recommended to also set
options(future.globals.maxSize= +Inf).

## Examples

``` r
#Create a VML data.frame
VMR_df <- data.frame(seqnames = c("chr21", "chr21"),
  start = c(10861376, 10862171),
  end = c(10862507, 10883548),
  probes =  I(list(c("cg15043638", "cg18287590", "cg17975851"),
                   c("cg13893907", "cg17035109", "cg06187584"))))

# Compute median correlation for each VMR
medCorVMR(VMR_df = VMR_df, methylation_data = RAMEN::test_methylation_data)
#>   seqnames    start      end       probes median_correlation
#> 1    chr21 10861376 10862507 cg150436....          0.7275160
#> 2    chr21 10862171 10883548 cg138939....          0.7468144


```
