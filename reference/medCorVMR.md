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
