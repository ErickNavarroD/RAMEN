# Identify Variable Methylated Regions in microarrays

Identifies autosomal Highly Variable Probes (HVP) and merges them into
Variable Methylated Regions (VMRs) given an Illumina manifest.

## Usage

``` r
findVMRs(
  array_manifest,
  methylation_data,
  cor_threshold = 0.15,
  var_method = "variance",
  var_threshold_percentile = 0.9,
  max_distance = 1000
)
```

## Arguments

- array_manifest:

  Information about the probes on the array. Requires the columns
  MAPINFO (basepair position of the probe in the genome), CHR
  (chromosome), TargetID (probe name) and STRAND (this is very important
  to set up, since the VMRs will only be created based on CpGs on the
  same strand; if the positions are reported based on a single DNA
  strand, this should contain either a vector of only "+", "-" or "\*"
  for all of the probes).

- methylation_data:

  A data frame containing M or B values, with samples as columns and
  probes as rows. Data is expected to have already passed through
  quality control and cleaning steps.

- cor_threshold:

  Numeric value (0-1) to be used as the median pearson correlation
  threshold for identifying VMRs (i.e. all VMRs will have a median
  pairwise probe correlation of this parameter).

- var_method:

  Method to use to measure variability in the data set. The options are
  "mad" (median absolute deviation) or "variance".

- var_threshold_percentile:

  The percentile (0-1) to be used as cutoff to define Highly Variable
  Probes (and therefore VMRs). The default is 0.9 because this
  percentile has been traditionally used in previous studies.

- max_distance:

  Maximum distance allowed for two probes to be grouped into a region.
  The default is 1000 because this window has been traditionally used in
  previous studies.

## Value

A list with the following elements:

- \$var_score_threshold: threshold used to define Highly Variable Probes
  (mad or variance, depending on the specified choice).

- \$highly_variable_probes: a data frame with the probes that passed the
  variability score threshold imposed by the user, and their variability
  score (MAD score or variance).

- \$canonical_VMRs: a GRanges object with strict candidate VMRs -
  regions composed of two or more contiguous, correlated and proximal
  Highly Variable Probes; thresholds depend on the ones specified by the
  user)

- \$non_canonical_VMRs: a GRanges object with highly variable probes
  without neighboring CpGs measured in *max_distance* on the array.
  Category created to take into acccount the Illumina array design of
  single probes capturing the methylation state of regulatory regions.

## Details

This function identifies HVPs using MAD scores or variance metrics, and
groups them into VMRs, which are defined as clusters of proximal and
correlated HVPs (distance and correlation defined by the user). Output
VMRs can be separated into canonical and non canonical. Canonical VMRs
are regions that meet the correlation and closeness criteria. For
guidance on which correlation threshold to use, we recommend checking
the Supplementary Figure 1 of the CoMeBack R package (Gatev *et al.*,
2020) where a simulation to empirically determine a default guidance
specification for a correlation threshold parameter dependent on sample
size is done. As default, we use a threshold of 0.15 as per the CoMeBack
authors minimum threshold suggestion. On the other hand, non canonical
VMRs are regions that are composed of HVPs that have no nearby probes
measured in the array (according to the max_distance parameter); this
category was created to account for the Illumina EPIC array design,
which has a high number of probes in regulatory regions that are
represented by a single probe. Furthermore, these probes have been shown
to be good representatives of the methylation state of its surroundings
(Pidsley et al., 2016). By creating this category, we recover those
informative HVPs that otherwise would be excluded from the analysis
because of the array design.

This function uses GenomicRanges::reduce() to group the regions, which
is strand-sensitive. In the Illumina microarrays, the MAPINFO for all
the probes is usually provided as for the + strand. If you are using
this array, we recommend to first convert the strand of all the probes
to "+".

This function supports parallel computing for increased speed. To do so,
you have to set the parallel backend in your R session BEFORE running
the function (e.g., doFuture::registerDoFuture()) and then the
evaluation strategy (e.g., future::plan(multisession)). After that, the
function can be run as usual. When working with big datasets, the
parallel backend might throw an error if you exceed the maximum allowed
size of globals exported for future expression. This can be fixed by
increasing the allowed size (e.g. running
options(future.globals.maxSize= +Inf) )

Note: this function excludes sex chromosomes.

## Examples

``` r
#We need to modify the RAMEN::test_array_manifest object by assigning to
#row names to the probe ID column; it was saved this way because storing
#the TargetID as row names reduced significantly the size of the data set.
test_array_manifest_final = RAMEN::test_array_manifest %>%
tibble::rownames_to_column(var = "TargetID")

VMRs = RAMEN::findVMRs(array_manifest = test_array_manifest_final,
                       methylation_data = RAMEN::test_methylation_data,
                       cor_threshold = 0,
                       var_method = "variance",
                       var_threshold_percentile = 0.9,
                       max_distance = 1000)
#> Identifying Highly Variable Probes...
#> Identifying non canonical Variable Methylated Regions...
#> Identifying canonical Variable Methylated Regions...
#> Applying correlation filter to canonical Variable Methylated Regions...
#> Warning: executing %dopar% sequentially: no parallel backend registered
```
