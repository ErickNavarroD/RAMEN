# Summarize the methylation states of Variable Methylated Loci (VML)

This function computes a representative methylation score for each
Variable Methylated Locus (VML) in a dataset. It returns a data frame
with the median methylation of each region per individual. For each VML
in a dataset, returns a with the median methylation of that region
(columns) per individual (rows) as representative score.

## Usage

``` r
summarizeVML(VML, methylation_data)
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

A matrix with samples as rows, and VML as columns. The value inside each
cell corresponds to the summarized methylation value of said VML in the
corresponding individual. The column names correspond to the VML_index.

## Details

This function supports parallel computing for increased speed. To do so,
you have to set the parallel backend in your R session BEFORE running
the function (e.g., *doParallel::registerDoParallel(4)*). After that,
the function can be run as usual.

## Examples

``` r
## Find VML in test data
# Set the parallel backend to use 2 workers
doParallel::registerDoParallel(2)
VML <- findVML(
  methylation_data = test_methylation_data,
  array_manifest = "IlluminaHumanMethylationEPICv1",
  cor_threshold = 0,
  var_method = "variance",
  var_distribution = "ultrastable",
  var_threshold_percentile = 0.99,
  max_distance = 1000
)
#> Identifying Highly Variable Probes...
#> Identifying sparse Variable Methylated Probes
#> Identifying Variable Methylated Regions...
#> Applying correlation filter to Variable Methylated Regions...

## Summarize methylation states of the found VML
summarized_VML <- summarizeVML(
  # Use only 5 for demonstration purposes
  VML = VML$VML[1:5, ],
  methylation_data = test_methylation_data
)
```
