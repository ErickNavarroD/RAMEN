# Summarize the methylation states of Variable Methylated Loci (VML)

This function computes a representative methylation score for each
Variable Methylated Locus (VML) in a dataset. It returns a data frame
with the median methylation of each region per individual. For each VML
in a dataset, returns a with the median methylation of that region
(columns) per individual (rows) as representative score.

## Usage

``` r
summarizeVML(VML_df, methylation_data)
```

## Arguments

- VML_df:

  A GRanges-like data frame. Must contain the following columns:
  "seqnames", "start", "end" and "probes" (containing lists as elements,
  where each contains a vector with the probes constituting the VML).
  This is the "VML" object returned by the *findVML()* function.

- methylation_data:

  A data frame containing M or B values, with samples as columns and
  probes as rows. Row names must be the CpG probe IDs.

## Value

A data frame with samples as rows, and VML as columns. The value inside
each cell corresponds to the summarized methylation value of said VML in
the corresponding individual. The column names correspond to the
VML_index.

## Details

This function supports parallel computing for increased speed. To do so,
you have to set the parallel backend in your R session BEFORE running
the function (e.g., *doParallel::registerDoParallel(4)*). After that,
the function can be run as usual.
