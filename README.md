
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RAMEN <a href="https://github.com/ErickNavarroD/RAMEN"><img src="man/figures/logo.png" align="right" height="150"/></a>

<!-- badges: start -->

[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/585986641.svg)](https://zenodo.org/badge/latestdoi/585986641)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Codecov test
coverage](https://codecov.io/gh/ErickNavarroD/RAMEN/graph/badge.svg)](https://app.codecov.io/gh/ErickNavarroD/RAMEN)
[![R-CMD-check](https://github.com/ErickNavarroD/RAMEN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ErickNavarroD/RAMEN/actions/workflows/R-CMD-check.yaml)
[![Status at rOpenSci Software Peer
Review](https://badges.ropensci.org/743_status.svg)](https://github.com/ropensci/software-review/issues/743)
<!-- badges: end -->

## Overview

Regional Association of Methylome variability with the Exposome and
geNome (RAMEN) is an R package which goal is to estimate the
contribution of genetic variants and environmental exposures to DNA
methylation variability at a genome-wide scale. RAMEN takes advantage of
the fact that DNA methylation levels at nearby CpG sites are often
correlated, and uses this information to identify Variable Methylated
Loci (VML) from microarray DNA methylation data. Then, using genomic and
exposomic data, it can identify which model out of the following
explains best the DNA methylation variability at each VML: genetic (G),
environmental (E), additive (G+E) or interactive (GxE).

RAMEN provides a Findable, Accesible, Interoperable and Reusable (FAIR)
workflow to conduct gene-environment contribution analyses to
high-dimensional DNA methylome data (described in [Navarro-Delgado et
al. (2025)](https://doi.org/10.1186/s13059-025-03864-4). Using a blend
of traditional statistical methods and machine learning approaches,
RAMEN is designed to be computationally efficient and user-friendly,
allowing researchers to gain insights into the complex interplay between
genetics, environment and DNA methylation variability. The package
includes a detailed tutorial, and individual functions that could be
useful for other applications beyond the gene-environment contribution
analysis.

## Installation

You can install the development version of RAMEN from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ErickNavarroD/RAMEN")
```

## Usage

For a detailed tutorial on how to use RAMEN, please check the package’s
vignette or
[website](https://ericknavarrod.github.io/RAMEN/articles/RAMEN.html).
Altogether, RAMEN provides a workflow that takes a set of individuals
with genome, exposome and DNA methylome information, and generates an
estimation of the contribution of genetic variants and environmental
exposures to its DNA methylation variability. Functions that conduct
computationally intensive tasks are compatible with parallel computing.

<img src="man/figures/RAMEN_pipeline.png" width="600"/>

In brief, the standard workflow consists of the following steps:

1.  Identify Variable Methylated Loci (VML) with `findVML()`.

``` r
library(RAMEN)
#>          __      _             ___
#>          )_)    /_)    )\/)    )_     )\  )
#>         / \    / /    (  (    (__    (  \(
#> 
#>                   (   )  (  (
#>                  (  ( )  (  )
#>                    )    )  (
#>                 _.(--'(''--.._
#>                /, _..-----).._,\
#>               |  `'''-----'''`  |
#>                \               /
#>                 '.           .'
#>                   '--.....--'
#> 
#> If you use RAMEN for your analysis, please cite: Navarro-Delgado, E.I.,
#> Czamara, D., Edwards, K. et al. RAMEN: Dissecting individual, additive
#> and interactive gene-environment contributions to DNA methylome variability
#> in cord blood. Genome Biol 26, 421 (2025).
#> https://doi.org/10.1186/s13059-025-03864-4
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)

VML <- RAMEN::findVML(
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

head(VML$VML) # Take a look at the identified VML data frame
#>   VML_index type seqnames    start      end width strand       probes n_VMPs
#> 1      VML1  VMR    chr21 10990119 10990903   785      + cg098720....      2
#> 2      VML2  VMR    chr21 11109021 11109336   316      + cg007508....      2
#> 3      VML3  VMR    chr21 31799091 31799248   158      + cg245007....      2
#> 4      VML4  VMR    chr21 32715908 32716792   885      + cg164170....      2
#> 5      VML5  VMR    chr21 15955548 15955699   152      - cg147721....      2
#> 6      VML6  VMR    chr21 26573136 26573196    61      - cg111120....      2
#>   median_correlation
#> 1          0.6099180
#> 2          0.6261681
#> 3          0.7279154
#> 4          0.6932442
#> 5          0.8120654
#> 6          0.6173683
```

2.  Summarize the regional methylation state of each VML with
    `summarizeVML()`.

``` r
summarized_methyl_VML <- RAMEN::summarizeVML(
  VML_df = VML$VML,
  methylation_data = test_methylation_data
)

# Look at the resulting object
summarized_methyl_VML[1:5, 1:5]
#>         VML1     VML2     VML3     VML4     VML5
#> ID1 4.935942 2.853168 6.389600 9.017997 2.714379
#> ID2 1.879166 2.699689 7.790474 3.134218 2.223942
#> ID3 3.311818 1.078262 4.135771 2.864724 8.648046
#> ID4 6.558106 4.683173 6.153156 3.828411 1.448140
#> ID5 2.899969 4.930614 4.919235 3.664651 2.926548
```

3.  Identify the SNPs in *cis* of each VML with `findCisSNPs()`.

``` r
VML_cis_snps <- RAMEN::findCisSNPs(
  VML_df = VML$VML,
  genotype_information = RAMEN::test_genotype_information,
  distance = 1e+06
)
#> Reminder: please make sure that the positions of the VML data frame and the ones in the genotype information are from the same genome build.

# Take a look at the result
head(VML_cis_snps)
#>   VML_index type seqnames    start      end width strand       probes n_VMPs
#> 1      VML1  VMR    chr21 10990119 10990903   785      + cg098720....      2
#> 2      VML2  VMR    chr21 11109021 11109336   316      + cg007508....      2
#> 3      VML3  VMR    chr21 31799091 31799248   158      + cg245007....      2
#> 4      VML4  VMR    chr21 32715908 32716792   885      + cg164170....      2
#> 5      VML5  VMR    chr21 15955548 15955699   152      - cg147721....      2
#> 6      VML6  VMR    chr21 26573136 26573196    61      - cg111120....      2
#>   median_correlation surrounding_SNPs          SNP
#> 1          0.6099180                1 21:10873....
#> 2          0.6261681                1 21:10873....
#> 3          0.7279154              659 21:30813....
#> 4          0.6932442              855 21:31718....
#> 5          0.8120654              726 21:14957....
#> 6          0.6173683              788 21:25582....
```

4.  Conduct a LASSO-based feature selection strategy to identify
    potentially relevant *cis* SNPs and environmental variables with
    `selectVariables()`.

``` r
selected_variables <- RAMEN::selectVariables(
  VML_df = VML_cis_snps,
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix = RAMEN::test_environmental_matrix,
  covariates = RAMEN::test_covariates,
  summarized_methyl_VML = summarized_methyl_VML,
  seed = 1
)
#> Loading required package: foreach
#> Loading required package: rngtools

head(selected_variables)
#>   VML_index selected_genot selected_env
#> 1      VML1   21:10873.... E43, E3,....
#> 2      VML2   21:10873....          E43
#> 3      VML3   21:32782.... E15, E25....
#> 4      VML4                         E43
#> 5      VML5   21:15248....     E40, E43
#> 6      VML6   21:25648....          E43
```

5.  Fit linear single-variable genetic (G), environmental (E), pairwise
    additive (G+E) and pairwise interaction (GxE) linear models, and
    select the best explanatory model for each VML with `lmGE()`.

``` r
lmge_res <- RAMEN::lmGE(
  selected_variables = selected_variables,
  summarized_methyl_VML = summarized_methyl_VML,
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix = RAMEN::test_environmental_matrix,
  covariates = RAMEN::test_covariates,
  model_selection = "AIC"
)

# Check the output
head(lmge_res)
#>   VML_index model_group    variables tot_r_squared g_r_squared e_r_squared
#> 1      VML1         G+E 21:10873....     0.5527507   0.1996156   0.3420620
#> 2      VML2         G+E 21:10873....     0.5115764   0.2095645   0.2714356
#> 3      VML3         GxE 21:32782....     0.5680092   0.2717880   0.2039656
#> 4      VML4           E          E43     0.3008107          NA   0.2282407
#> 5      VML5         G+E 21:15671....     0.7548714   0.4246781   0.2232954
#> 6      VML6         G+E 21:27306....     0.5885708   0.2317301   0.1985512
#>   gxe_r_squared      AIC second_winner delta_aic delta_r_squared basal_AIC
#> 1            NA 144.6841           GxE 1.0497027    -0.013945293  164.4893
#> 2            NA 148.7250           GxE 1.3455746    -0.010539189  165.2905
#> 3    0.06543422 148.7964           G+E 0.9072693     0.043959423  167.1613
#> 4            NA 156.8403          <NA>        NA              NA  163.3152
#> 5            NA 131.9394           GxE 1.5262724    -0.003840405  166.7269
#> 6            NA 141.4738           GxE 1.6975735    -0.004126734  158.9477
#>   basal_rsquared
#> 1     0.01107306
#> 2     0.03057637
#> 3     0.02682137
#> 4     0.07256996
#> 5     0.10689789
#> 6     0.15828947
```

6.  Simulate a null distribution of G and E effects on DNAme variability
    with `nullDistGE()`, and use it to filter out poor-performing best
    explanatory models selected by *lmGE()*.

``` r
null_dist <- RAMEN::nullDistGE(
  VML_df = VML_cis_snps,
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix = RAMEN::test_environmental_matrix,
  summarized_methyl_VML = summarized_methyl_VML,
  permutations = 1,
  covariates = RAMEN::test_covariates,
  seed = 1,
  model_selection = "AIC"
)
#> Starting permutation 1 of 1
#> Starting variable selection of permutation 1 of 1
#> Starting lmGE in permutation 1 of 1
#> Wrapping up permutation 1 of 1

#Set threshold
cutoff_single <- quantile(
  null_dist %>%
    filter(model_group %in% c("G", "E")) %>%
    pull(R2_difference),
  0.95
)
cutoff_joint <- quantile(
  null_dist %>%
    filter(model_group %in% c("G+E", "GxE")) %>%
    pull(R2_difference),
  0.95
)

# Get a data frame with the final results
final_res <- lmge_res %>%
  dplyr::mutate(
    r2_difference_basal = tot_r_squared - basal_rsquared,
    # Label if the best explanatory model passes its corresponding threshold
    pass_cutoff_threshold = case_when(
      model_group %in% c("G", "E") ~ r2_difference_basal > cutoff_single,
      model_group %in% c("G+E", "GxE") ~ r2_difference_basal > cutoff_joint
    ),
    # Label the final model group, replacing bad performing winning models with
    #  "B" (basal)
    model_group = case_when(
      pass_cutoff_threshold ~ model_group,
      TRUE ~ "B"
    )
  ) %>%
  dplyr::select(-pass_cutoff_threshold) # Drop temporary column

# Keep only VML that have informative models with out data
filtered_res <- final_res %>%
  dplyr::filter(!model_group == "B") # Filter based on the cutoff threshold

final_res %>%
  dplyr::group_by(model_group) %>%
  dplyr::summarise(count = n()) %>%
  ggplot2::ggplot(aes(x = model_group, y = count)) +
  ggplot2::geom_col() +
  ggplot2::xlab("Best explanatory model") +
  ggplot2::ylab("VML") +
  ggplot2::theme_classic()
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

## Variations to the standard workflow

Besides using RAMEN for a gene-environment contribution analysis, the
package provides individual functions that could help users in other
tasks, such as:

- Reduction of multiple hypothesis test burden in EWAS or differential
  methylation analysis by using VML instead of individual probes.
- Fit additive and interaction models given a set of variables of
  interest and select the best explanatory model for DNAme data
  (e.g. epistasis or ExE studies).
- Quickly identify SNPs in *cis* of CpG probes.
- Get the median correlation of probes in custom regions of interest
  with `medCorVMR()`.

## How to get help for RAMEN

If you have any question about RAMEN usage, please [post a new
issue](https://github.com/ErickNavarroD/RAMEN/issues/new/choose) in this
github repository so that future users also benefit from the discussion.

## Acknowledgments

This package was developed by Erick I. Navarro-Delgado under the
supervision of Dr. Keegan Korthauer and Dr. Michael S. Kobor. We want to
thank the members of the Kobor and Korthauer lab for their feeback
during the development of RAMEN. Additionally, we want to thank Carlos
Cortés-Quiñones and Dorothy Lin for helping create the package logo.
Erick conceptualized the logo, Carlos drew it, and Dorothy refined it
and finished the lettering.

## Funding

This work was supported by the University of British Columbia, the BC
Children’s Hospital Research Institute and the Social Exposome Cluster.

## Citing RAMEN

If you use RAMEN for any of your analyses, please cite the following
publication:

- Navarro-Delgado, E.I., Czamara, D., Edwards, K. et al. RAMEN:
  Dissecting individual, additive and interactive gene-environment
  contributions to DNA methylome variability in cord blood. *Genome
  Biol* 26, 421 (2025). <https://doi.org/10.1186/s13059-025-03864-4>

## Code of conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

## Licence

GPL (\>= 3)
