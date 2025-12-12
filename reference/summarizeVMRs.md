# Summarize the methylation states of Variable Methylated Regions (VMRs)

For each VMR in a dataset, returns an object with the median methylation
of that region per individual as representative score.

## Usage

``` r
summarizeVMRs(VMRs_df, methylation_data)
```

## Arguments

- VMRs_df:

  A GRanges object converted to a data frame. Must contain the following
  columns: "seqnames", "start", "end" (all of which are produced
  automatically when doing the object conversion) and "probes"
  (containing a list where each element contains a vector with the
  probes constituting the VMR).

- methylation_data:

  A data frame containing M or B values, with samples as columns and
  probes as rows.

## Value

A data frame with samples as rows, and VMRs as columns. The value inside
each cell corresponds to the summarized methylation value of said VMR in
the corresponding individual. The column names correspond to the
VMR_index, which is created if not already existing based on the
rownames of the VMR_df.

## Details

This function supports parallel computing for increased speed. To do so,
you have to set the parallel backend in your R session BEFORE running
the function (e.g., doFuture::registerDoFuture()) and then the
evaluation strategy (e.g., future::plan(multisession)). After that, the
function can be run as usual.
