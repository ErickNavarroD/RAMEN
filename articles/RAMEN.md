# RAMEN

## Introduction

**Regional Association of Methylome variability with the Exposome and
geNome (RAMEN)** is an R package whose goal is to integrate genomic,
methylomic and exposomic data to model the contribution of genetics (G)
and the environment (E) to DNA methylation (DNAme) variability. RAMEN
identifies Variable Methylated Regions (VMRs) in microarray DNAme data
and then, using genotype and environmental data, it identifies which of
the following models better explains this variability in regions across
the methylome:

| Model                         | Name                   | Abbreviation |
|:------------------------------|:-----------------------|:-------------|
| DNAme ~ G + covars            | Genetics               | G            |
| DNAme ~ E + covars            | Environmental exposure | E            |
| DNAme ~ G + E + covars        | Additive               | G+E          |
| DNAme ~ G + E + G\*E + covars | Interaction            | GxE          |

Fitted models

where G variables are represented by SNPs, E variables by environmental
exposures, and where covars are concomitant variables (i.e. variables
that are adjusted for in the model and not of interest in the study such
as cell type proportion, age, etc.).

The main [gene-environment interaction
modeling](#gene-environment-interaction-analysis) pipeline is conducted
though six core functions:

- [`findVMRs()`](https://ericknavarrod.github.io/RAMEN/reference/findVMRs.md)
  identifies Variable Methylated Regions (VMRs) in microarrays
- [`summarizeVMRs()`](https://ericknavarrod.github.io/RAMEN/reference/summarizeVMRs.md)summarizes
  the regional methylation state of each VMR
- [`findCisSNPs()`](https://ericknavarrod.github.io/RAMEN/reference/findCisSNPs.md)
  identifies the SNPs in *cis* of each VMR
- [`selectVariables()`](https://ericknavarrod.github.io/RAMEN/reference/selectVariables.md)
  conducts a LASSO-based variable selection strategy to identify
  potentially relevant *cis* SNPs and environmental variables
- [`lmGE()`](https://ericknavarrod.github.io/RAMEN/reference/lmGE.md)
  fits linear single-variable genetic (G) and environmental (E), and
  pairwise additive (G+E) and interaction (GxE) linear models and select
  the best explanatory model per VMR.
- [`nullDistGE()`](https://ericknavarrod.github.io/RAMEN/reference/nullDistGE.md)
  simulates a delta R squared null distribution of G and E effects on
  DNAme variability. Useful for filtering out poor-performing best
  explanatory models selected by *lmGE()*.

These functions are compatible with parallel computing, which is
recommended due to the computationally intensive tasks conducted by the
package.

In addition to the [standard gene-environment interaction modeling
pipeline](#gene-environment-interaction-analysis), RAMEN can be useful
for other DNAme analyses (see [Variations to the standard
workflow](#variations-to-the-standard-workflow)), such as reducing the
tested sites in Epigenome Wide Association Studies, grouping DNAme
probes into regions, identifying SNPs near a probe, etc.

### Citation

The manuscript detailing RAMEN and its use is currently under
preparation. For more information about this please contact Erick I.
Navarro-Delgado at <erick.navarrodelgado@bcchr.ca>.

## Gene-environment interaction analysis

The main purpose of the RAMEN package is to conduct a methylome-wide
analysis to identify which model (G, E, G+E or GxE) better explains the
variability across the genome. In this vignette, we will illustrate how
to use the package.

To conduct this analysis, the following cleaned data sets (i.e. after
quality and exploratory data analysis checks) from a cohort are
required:

- DNAme data
- DNAme array manifest
- Genotyping data
- Genotype information
- Environmental exposure data
- Concomitant variables data

Once that we have that data, the overview of the pipeline is the
following:

![RAMEN pipeline](RAMEN_pipeline.png)

RAMEN pipeline

where:

- DNAme data is grouped into VMRs, and then the DNAme state per
  individual is summarized in each region
- Using the identified VMRs and the genomic information, we identify the
  SNPs in *cis* for each VMR
- Both the *cis* SNPs and the exposome data are subjected to the
  variable selection stage
- The selected variables (SNPs a and Es) enter the modelling stage,
  which outputs one single winning model per VMR
- The thresholds obtained from the simulated null distribution are used
  to remove winning models which performance are likely to be due to
  chance.

In the following sections we will go through each of these steps and
guide the user regarding the recommended parameters to use in each
function of the package. For illustration purposes, we provide small toy
data sets that do not intend to simulate a real biological phenomenon.
These data sets are already available in the RAMEN package.

``` r
#Load the packages used throughout the vignette
library(RAMEN)
library(dplyr)
library(ggplot2)
library(tidyr)
```

### Identify VMRs and summarize their methylation state

The first step of the pipeline is to identify the **Variable Methylated
Regions**(VMRs) in the data set. You might be wondering *“What is a VMR
and why do we use them instead of DNAme levels from each CpG site?”*. We
use **regions** because it is well established that nearby CpG sites are
[very likely to share a similar DNAme
profile](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6093082/) and
therefore work as functional units. Then, from a statistical point of
view, testing separately proximal CpGs that are part of the same unit is
redundant. On the other side, we use only **variable** regions because
we are interested in the units that display a high level of variability;
in other words, in non-variant sites there is no variability left to be
explained by genetics or environment. So, in conclusion, we use **VMRs**
to increase our power and reduce the multiple hypothesis testing burden
by grouping probes that are likely to work as a biological unit, and by
only focusing in the set of regions that are of interest of this study.

RAMEN identifies 2 categories of VMRs:

- Canonical VMRs: Group of Highly Variable Probes that are proximal and
  correlated. Highly Variable Probes are defined as probes above a
  specific variance percentile threshold specified by the user (the
  default is 90th percentile). The proximity distance and pearson
  correlation threshold is specified by the user, and the defaults are 1
  kilobase and 0.15 respectively. For guidance on which correlation
  threshold to use, we recommend checking the Supplementary Figure 1 of
  the CoMeBack R package (Gatev et al., 2020) where a simulation to
  empirically determine a default guidance specification for a
  correlation threshold parameter dependent on sample size is done.
- Non canonical VMRs: Regions that are composed of Highly Variable
  Probes that have no nearby probes measured in the array (according to
  the distance parameter specified by the user). This category was
  created to take into account the characteristics of the DNAme
  microarray plataform, which has probes covering non-homogenelously the
  genome. This is specially important for microarrays such as the EPIC
  array which has a high number of probes in regulatory regions that are
  represented by a single probe. Furthermore, these probes have been
  shown to be good representatives of the methylation state of its
  surroundings (Pidsley et al., 2016). By creating this category, we
  recover those informative HVPs that would otherwise be excluded from
  the analysis because of working with the canonical VMR definition in
  the array context.

The first step is to identify **Variable Methylated Regions**(VMRs)
using the
[`RAMEN::findVMRs()`](https://ericknavarrod.github.io/RAMEN/reference/findVMRs.md)
function. This function uses GenomicRanges::reduce() to group the
regions, which is strand-sensitive. In the Illumina microarrays, the
MAPINFO for all the probes is usually provided as for the + strand. If
you are using this array, we recommend to first convert the strand of
all the probes to “+”.

For this step, we recommend users to use M-values because its use is
more appropiate for statistical analyses (see Pan Du, *et al.*, 2010,
*BMC Bioinformatics*)

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

#Take a look at the resulting object 
dplyr::glimpse(VMRs)
#> List of 4
#>  $ var_score_threshold   : Named num 13.5
#>   ..- attr(*, "names")= chr "90%"
#>  $ highly_variable_probes:'data.frame':  300 obs. of  2 variables:
#>   ..$ TargetID : chr [1:300] "cg06187584" "cg09872009" "cg05437132" "cg26689799" ...
#>   ..$ var_score: num [1:300] 17.1 14.8 15.2 13.7 17.1 ...
#>  $ canonical_VMRs        :Formal class 'GRanges' [package "GenomicRanges"] with 7 slots
#>   .. ..@ seqnames       :Formal class 'Rle' [package "S4Vectors"] with 4 slots
#>   .. ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
#>   .. ..@ strand         :Formal class 'Rle' [package "S4Vectors"] with 4 slots
#>   .. ..@ seqinfo        :Formal class 'Seqinfo' [package "Seqinfo"] with 4 slots
#>   .. ..@ elementMetadata:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
#>   .. ..@ elementType    : chr "ANY"
#>   .. ..@ metadata       : list()
#>  $ non_canonical_VMRs    :Formal class 'GRanges' [package "GenomicRanges"] with 7 slots
#>   .. ..@ seqnames       :Formal class 'Rle' [package "S4Vectors"] with 4 slots
#>   .. ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
#>   .. ..@ strand         :Formal class 'Rle' [package "S4Vectors"] with 4 slots
#>   .. ..@ seqinfo        :Formal class 'Seqinfo' [package "Seqinfo"] with 4 slots
#>   .. ..@ elementMetadata:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
#>   .. ..@ elementType    : chr "ANY"
#>   .. ..@ metadata       : list()
```

As we can see,
[`RAMEN::findVMRs()`](https://ericknavarrod.github.io/RAMEN/reference/findVMRs.md)
returns a list with four elements:

- var_score_threshold: The MAD-score or variance threshold used to
  define Highly Variable Probes.
- highly_variable_probes: a data frame with the probes that passed the
  variability score threshold imposed by the user, and their variability
  score (MAD score or variance).
- canonical_VMRs: a GRanges object with strict candidate VMRs - regions
  composed of two or more contiguous, correlated and proximal Highly
  Variable Probes; thresholds depend on the ones specified by the user)
- non_canonical_VMRs: a GRanges object with highly variable probes
  without neighboring CpGs measured in max_distance on the array.
  Category created to take into acccount the Illumina array design of
  single probes capturing the methylation state of regulatory regions.

Furthermore, we can see the following warning message in the chunk
above:

``` r
#> Warning: executing %dopar% sequentially: no parallel backend registered
```

This is printed in the screen just to warn us that
[`RAMEN::findVMRs()`](https://ericknavarrod.github.io/RAMEN/reference/findVMRs.md)
is running sequentially. RAMEN supports parallel computing for increased
speed, which is really important when working with real data sets that
tend to contain information from the whole genome. To do so, you have to
set the parallel backend in your R session BEFORE running the function
(e.g., `doFuture::registerDoFuture()`) and then the evaluation strategy
(e.g., `future::plan(multisession)`). After that, the function can be
run normally. When working with big datasets, the parallel backend might
throw an error if you exceed the maximum allowed size of globals
exported for future expression. This can be fixed by increasing the
allowed size (e.g. running `options(future.globals.maxSize= +Inf)`)

Finally, we will convert the output of
[`RAMEN::findVMRs()`](https://ericknavarrod.github.io/RAMEN/reference/findVMRs.md)
to a data frame, which is an object that we can easily use to produce
plots and explore the results, and the object that is needed for the
following parts of the pipeline.

``` r
VMRs_df = data.frame(VMRs[["canonical_VMRs"]]) %>% 
  rbind(as.data.frame(VMRs[["non_canonical_VMRs"]])) %>% 
  dplyr::select( -c(width.1,strand))

head(VMRs_df)
#>   seqnames    start      end width n_VMPs       probes median_correlation
#> 1       21 10990119 10990903   785      2 cg098720....          0.6099180
#> 2       21 11109021 11109336   316      2 cg007508....          0.6261681
#> 3       21 15456605 15456789   185      2 cg195094....          0.8595934
#> 4       21 15588038 15588830   793      2 cg023649....          0.6829510
#> 5       21 15955548 15956048   501      3 cg147721....          0.7548517
#> 6       21 16436438 16437768  1331      3 cg116304....          0.7124538
```

With the VMRs as a data frame, we can explore our data set using ggplot,
such as the following example:

``` r
VMRs_df %>% 
  dplyr::filter(width > 1) %>% #Only plot canonical VMRs 
  ggplot2::ggplot(aes(x = width))+
  ggplot2::geom_histogram(binwidth = 50, fill = "#BAB4D8")+
  ggplot2::theme_classic()+
  ggplot2::ggtitle("Canonical VMRs width (bp)") 
```

![Width of non canonical VMRs (base
pairs).](RAMEN_files/figure-html/unnamed-chunk-6-1.png)

Width of non canonical VMRs (base pairs).

Next, we want to summarize the DNAme level of each VMR per individual.
To do this, we use
[`RAMEN::summarizeVMRs()`](https://ericknavarrod.github.io/RAMEN/reference/summarizeVMRs.md).
For non canonical VMRs, there is nothing to summarize, so the DNAme
level of the corresponding probe is returned. For canonical VMRs, the
median DNAme level of the region is returned per individual as the
representative value.

``` r
summarized_methyl_VMR = RAMEN::summarizeVMRs(VMRs_df = VMRs_df,
                                             methylation_data = test_methylation_data)

# Look at the resulting object
summarized_methyl_VMR[1:5,1:5]
#>            1        2        3        4           5
#> ID1 4.935942 2.853168 3.219038 5.384813 3.320883317
#> ID2 1.879166 2.699689 3.952487 2.927166 1.934283315
#> ID3 3.311818 1.078262 1.153573 4.690962 9.238518404
#> ID4 6.558106 4.683173 1.688233 3.448903 0.007576317
#> ID5 2.899969 4.930614 2.243714 5.788151 2.843904388
```

The result is a data frame of VMR IDs as columns and individual IDs as
rows.

### Identify *cis* SNPs

After identifying the VMRs, we recommend to use only SNPs in *cis* of
each VMR, since genetic variants that associate with DNAme changes tend
to be more abundant in the surroundings of the corresponding DNAme site
(McClay *et al.*, 2015). Also, the effect sizes of mQTLs (genetic
variants associated with DNAme changes) are stronger in *cis* SNPs
compared to *trans* SNPs. Then, by restricting the analysis to *cis*
SNPs, we greatly reduce the number of variables while keeping most of
the important ones.

There is not a clear consensus on how close a SNP has to be from a DNAme
site to be considered *cis* - the distance threshold tend to go from few
kb to 1 megabase. We recommend to use a 1 Mb window to cast a wide net
to catch potentially relevant SNPs.

``` r
VMRs_df = RAMEN::findCisSNPs(VMRs_df = VMRs_df, 
                            genotype_information = RAMEN::test_genotype_information, 
                            distance = 1e+06)
#> Important: please make sure that the positions of the VMR data frame and the ones in the genotype information are from the same genome build.

#Take a look at the result
dplyr::glimpse(VMRs_df) 
#> Rows: 131
#> Columns: 10
#> $ VMR_index          <chr> "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", …
#> $ seqnames           <fct> 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,…
#> $ start              <int> 10990119, 11109021, 15456605, 15588038, 15955548, 1…
#> $ end                <int> 10990903, 11109336, 15456789, 15588830, 15956048, 1…
#> $ width              <int> 785, 316, 185, 793, 501, 1331, 888, 977, 573, 61, 5…
#> $ n_VMPs             <int> 2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, …
#> $ probes             <list> <"cg09872009", "cg05437132">, <"cg00750806", "cg12…
#> $ median_correlation <dbl> 0.6099180, 0.6261681, 0.8595934, 0.6829510, 0.75485…
#> $ surrounding_SNPs   <int> 1, 1, 610, 667, 726, 794, 757, 529, 709, 788, 829, …
#> $ SNP                <list> "21:10873592:G:A", "21:10873592:G:A", <"21:1459574…
```

We can see that the resulting VMRs_df object is almost exactly the same,
but with two new columns (*surrounding_SNPs* and *SNP*) that contain
information about how many SNPs were found in *cis* and what are their
IDs according to the genotype data that we have.

It is important to highlight the following characteristics of the
VMRs_df object:

- The *VMR_index* column is a character vector. This column corresponds
  to the unique identifier of each VMR in our data set. It is important
  to **keep it as a character** and be careful with it not being
  converted to numeric, which can happen if you save the VMRs_df object
  as a table, read it, and use that second object in the rest of the
  pipeline.
- The columns *probes* and *SNP* contain **lists**. This structure is
  really important for the rest of the analysis and columns containing
  lists will keep appearing in other function outputs. If you want to
  know the recommended way to save and load these objects, please check
  the [Frequently Asked Questions](#frequently-asked-questions).

We can also explore the resulting object through plots such as the
following:

``` r
VMRs_df %>% 
  dplyr::mutate(type = case_when(n_VMPs == 1 ~ "non canonical", #non canonical VMRs have only 1 probe by definition
                          TRUE ~ "canonical"),
         surrounding_SNPs = case_when( surrounding_SNPs > 3000 ~ 3000,
                                       TRUE ~ surrounding_SNPs)) %>% 
  ggplot2::ggplot(aes(x = surrounding_SNPs)) +
  ggplot2::geom_density() +
  ggplot2::facet_grid("type") +
  ggplot2::theme_classic()
```

![Disribution of SNPs in cis of
VMRs.](RAMEN_files/figure-html/cissnps-1.png)

Disribution of SNPs in cis of VMRs.

### Conduct variable selection on genome and exposome variables

The following stage in the pipeline is to screen the available variables
in our environmental and *cis* SNPs data sets to identify the
potentially relevant ones. This is achieved with the
[`RAMEN::selectVariables()`](https://ericknavarrod.github.io/RAMEN/reference/selectVariables.md)
function. This function uses a data-driven approach based on LASSO,
which is an embedded variable selection method commonly used in machine
learning.

In a nutshell, LASSO penalizes models that are more complex (i.e., that
contain more variables) in favor of simpler models (i.e. that contain
less variables), but not at the expense of reducing predictive power.
Using LASSO’s variable screening property (i.e., with high probability,
the LASSO estimated model includes the substantial covariates and drops
the redundant ones) this function is intended to select genotype and
environmental variables in each Variable Methylated Region (VMR) with
potential relevance in the presence of the user-specified concomitant
variables (which are known DNAme confounders such as age, cell type
proportion, etc.). For more information about the method, we encourage
the users to read the documentation of the function, and for further
information about LASSO we direct readers to Bühlmann and Van de Geer,
2011.

Overall, conducting our variable selection strategy reduces the overall
computational time and improves the modeling performance by:

- Reducing the universe of variables that will be used to fit models in
  the following stage (G/E/G+E/GxE model fitting and comparison)
- Removing redundant variables, which are highly expected in genetic and
  environmental data sets with a high number of variables
- Limiting the interactions terms to scenarios where both the G and E
  main effects were selected to be potentially relevant, which can be
  think of as an interaction variable selection using a weak hierarchy
- Using LASSO, a method with good variable selection performance and
  scalability

Please make sure that your data has no NAs, since the LASSO
implementation we use in RAMEN does not support missing values. If your
data has missing values, consider
[handling](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3668100/) them.

``` r
selected_variables = RAMEN::selectVariables(
  VMRs_df = VMRs_df,
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix= RAMEN::test_environmental_matrix,
  covariates = RAMEN::test_covariates,
  summarized_methyl_VMR = summarized_methyl_VMR,
  seed = 1
)
#> Loading required package: foreach
#> Loading required package: rngtools
```

Since LASSO makes use of Random Number Generation, setting a seed is
highly encouraged for result’s reproducibility using the *seed*
argument. As a note, setting a seed inside of this function modifies the
seed globally (which is R’s default behavior).

The output of
[`RAMEN::selectVariables()`](https://ericknavarrod.github.io/RAMEN/reference/selectVariables.md)
is an object with the VMR index, and the G and E variables selected for
each VMR.

``` r
dplyr::glimpse(selected_variables)
#> Rows: 131
#> Columns: 3
#> $ VMR_index      <chr> "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"…
#> $ selected_genot <list> "21:10873592:G:A", "21:10873592:G:A", <"21:15450190:C:…
#> $ selected_env   <list> <"E43", "E3", "E5", "E7", "E24", "E25", "E28", "E35", …
```

We can see how using
[`RAMEN::selectVariables()`](https://ericknavarrod.github.io/RAMEN/reference/selectVariables.md)
reduces the number of variables (originally 100 environmental variables
and 785 SNPs per VMR on average as seen in Figure @ref(fig:cissnps)).

``` r
selected_variables %>% 
  dplyr::left_join(VMRs_df %>% 
              select(c(VMR_index,n_VMPs)),
            by = "VMR_index") %>% 
  dplyr::transmute(VMR_index = VMR_index,
            VMR_type = case_when(n_VMPs > 1 ~ "canonical",
                          n_VMPs == 1 ~ "non canonical"), 
            genotype = lengths(selected_genot),
            environment = lengths(selected_env)) %>% 
  tidyr::pivot_longer(-c(VMR_index, VMR_type)) %>% 
  dplyr::rename(group = name,
         variables = value) %>% 
  ggplot2::ggplot(aes(x = VMR_type, y = variables)) +
  ggplot2::geom_violin() + 
  ggplot2::geom_boxplot(width=0.1, outlier.shape=NA) +
  ggplot2::facet_wrap(~group)+
  ggplot2::ggtitle("Selected variables")  +
  ggplot2::theme_classic()
```

![Number of G and E selected
variables.](RAMEN_files/figure-html/selectedvars-1.png)

Number of G and E selected variables.

It is also expected in real data to have VMRs where no SNP and/or no
environmental variables were selected, since not all the DNAme sites in
the genome are expected to show an association with the genetic
variation or environmental exposures data sets that are captured in a
study. The proportion of VMRs under these scenarios will depend on the
data sets.

#### Author’s note about variables interpretation

LASSO variable selection is not consistent when there is
multicollinearity in the data (i.e., correlation between variables),
which is expected due to the high amount of G and E variables that are
present in studies of this kind. This means that if you were to run
LASSO several times, and two variables were to be highly correlated, the
method would select one and drop the other one at random. This is not a
problem with the pipeline because the main conclusion per VMR is whether
the DNAme is better explained by G and/or E components. As an example,
if a VMR is better explained by SNP1 and SNP2, which are both highly
correlated one with the other, LASSO will randomly pick SNP1 OR SNP2
(because they are relevant but they provide redundant information); if
we were to fit a model with SNP1 or SNP2 in the following stage, the
winning model would still be G. In other words, the main goal of the
pipeline is to know whether the VMR’s DNAme is better explained by G
and/or E. The user is therefore warned to **be cautious not to
over-interpret the individual selected variables**. Selected variables
might be used as hypothesis generators of associations, keeping in mind
that the selected variable might be representing other variables in the
data set that provide similar information.

### Identify the best explanatory model (G/E/G+E/GxE) per VMR

Now that we have selected the list of potentially relevant G and E
variables, we will fit the models mentioned in Table
@ref(tab:modelstable) using the
[`RAMEN::lmGE()`](https://ericknavarrod.github.io/RAMEN/reference/lmGE.md)
function. This function fits, for each VMR, G and E models with all of
the variables selected, as well as all their possible pairwise
combinations of G+E and GxE.

After fitting this model, the best model per group (group = G, E, G+E or
GxE) is selected using Akaike Information Criterion (AIC) or Bayesian
Information Criterion (BIC). We recommend using AIC because BIC assumes
that the true model is in the set of compared models. Since this
function fits models with individual variables, and we assume that DNAme
variability is more likely to be influenced by more than one single
SNP/environmental exposure at a time, we hypothesize that in most cases,
the true model will not be in the set of compared models. Also, AIC
excels in situations where all models in the model space are
“incorrect”, and AIC is preferentially used in cases where the true
underlying function is unknown and our selected model could belong to a
very large class of functions where the relationship could be pretty
complex. It is worth mentioning however that, both metrics tend to pick
the same model in a large number of scenarios. We suggest the users to
read Arijit Chakrabarti & Jayanta K. Ghosh, 2011 for further information
about the difference between these metrics. After selecting the best
model per group (G,E,G+E pr GxE), the model with the lowest AIC or BIC
will be declared as the winning model.

Additionally,
[`RAMEN::lmGE()`](https://ericknavarrod.github.io/RAMEN/reference/lmGE.md)
conducts a variance decomposition analysis, so that the relative R2
contribution of each of the variables of interest (G, E and GxE) is
reported. This decomposition is done using the
*[relaimpo](https://CRAN.R-project.org/package=relaimpo)* R package,
using the Lindeman, Merenda and Gold (lmg) method, which is based on the
heuristic approach of averaging the relative R contribution of each
variable over all input orders in the linear model.

``` r
lmge_res = RAMEN::lmGE(
  selected_variables = selected_variables,
  summarized_methyl_VMR = summarized_methyl_VMR,
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix = RAMEN::test_environmental_matrix,
  covariates = RAMEN::test_covariates,
  model_selection = "AIC"
)

# Check the output
dplyr::glimpse(lmge_res)
#> Rows: 128
#> Columns: 13
#> $ VMR_index       <chr> "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11…
#> $ model_group     <chr> "G+E", "G+E", "GxE", "GxE", "G+E", "GxE", "G+E", "GxE"…
#> $ variables       <list> <"21:10873592:G:A", "E43">, <"21:10873592:G:A", "E43"…
#> $ tot_r_squared   <dbl> 0.5527507, 0.5115764, 0.6075657, 0.6794845, 0.7680461,…
#> $ g_r_squared     <dbl> 0.1996156, 0.2095645, 0.3064024, 0.2738786, 0.3981840,…
#> $ e_r_squared     <dbl> 0.34206198, 0.27143559, 0.15293109, 0.20266291, 0.2769…
#> $ gxe_r_squared   <dbl> NA, NA, 0.1330680, 0.1126333, NA, 0.1432113, NA, 0.108…
#> $ AIC             <dbl> 144.6841, 148.7250, 145.6987, 135.7819, 129.9252, 123.…
#> $ second_winner   <chr> "GxE", "GxE", "G+E", "G+E", "GxE", "G+E", "GxE", "G+E"…
#> $ delta_aic       <dbl> 1.04970268, 1.34557464, 1.01022314, 3.88212651, 0.0761…
#> $ delta_r_squared <dbl> -0.013945293, -0.010539189, 0.041420496, 0.069427874, …
#> $ basal_AIC       <dbl> 164.4893, 165.2905, 167.3019, 161.0771, 166.8374, 159.…
#> $ basal_rsquared  <dbl> 0.011073059, 0.030576374, 0.015164242, 0.090309667, 0.…
```

The output of
[`RAMEN::lmGE()`](https://ericknavarrod.github.io/RAMEN/reference/lmGE.md)
is a data frame with the following 13 columns:

- *VMR_index*: The index of the respective VMR
- *model_group*: The selected winning model (G, E, G+E or GxE)
- *variables*: The variable(s) that are present in the winning model
  (excluding the covariates, which are included in all the models)
- *tot_r_squared*: total R squared of the winning model
- *g_r_squared*: Estimated R2 allocated to the G component in the
  winning model, if applicable
- *e_r_squared*: Estimated R2 allocated to the E in the winning model,
  if applicable.
- *gxe_r_squared*: Estimated R2 allocated to the interaction in the
  winning model (GxE), if applicable.
- *AIC/BIC*: AIC or BIC metric from the best model in each VMR
  (depending on the option specified in the argument model_selection).
- *second_winner*: The second group that possesses the next best model
  after the winning one (i.e., G, E, G+E or GxE). This column may have
  NA if the variables in selected_variables correspond only to one group
  (G or E), so that there is no other model groups to compare to.
- *delta_aic/delta_bic*: The difference of AIC or BIC (depending on the
  option specified in the argument model_selection) of the winning model
  and the best model from the second_winner group (i.e., G, E, G+E or
  GxE). This column may have NA if the variables in selected_variables
  correspond only to one group (G or E), so that there is no other
  groups to compare to.
- *delta_r_squared*: The R2 of the winning model - R2 of the second
  winner model. This column may have NA if the variables in
  selected_variables correspond only to one group (G or E), so that
  there is no other groups to compare to.
- *basal_AIC/basal_BIC*: AIC or BIC of the basal model (i.e., model
  fitted only with the concomitant variables specified in the covariates
  argument)
- *basal_rsquared*: The R2 of the basal model (i.e., model fitted only
  with the concomitant variables specified in the covariates argument)

#### Remove poor performing winning models

The core pipeline from the RAMEN package identifies the best explanatory
model per VMR. However, despite these models being winners in comparison
to models including any other G/E variable(s) in the dataset, some
winning models might perform no better than what we would expect by
chance. Therefore, The last step of the pipeline after identifying the
winning models per VMR is to compute a null distribution to remove the
winning models that are very likely to be winners by chance. To do so,
we use
[`RAMEN::nullDistGE()`](https://ericknavarrod.github.io/RAMEN/reference/nullDistGE.md).

The goal of
[`RAMEN::nullDistGE()`](https://ericknavarrod.github.io/RAMEN/reference/nullDistGE.md)
is to create a distribution of increase in R2 after including the
SNP/E/SNP\*E variables under the null hypothesis of G and E having no
associations with DNAme. The null distribution is obtained through
shuffling the G and E variables in a given dataset and conducting the
variable selection and G/E model selection. That way, we can simulate
how much additional variance would be explained by the models defined as
winners by the RAMEN methodology in a scenario where the G and E
associations with DNAme are randomized. This distribution can be then
used to filter out winning models in the non-shuffled dataset that do
not add more to the explained variance of the basal model than what
randomized data do.

For clarification, please note that in this vignette when we refer to
SNP\*E, we are referring to the interaction term that is present in the
the interaction model (i.e. interaction variable in the GxE model).

Under the assumption that after adjusting for the concomitant variables
all VMRs across the genome follow the same behavior regarding an
increment of explained variance with randomized G and E data, we can
pool the delta R squared values from all VMRs to create a null
distribution taking advantage of the high number of VMRs in the dataset.
This assumption decreases significantly the number of permutations
required to create a null distribution and reduces the computational
time. For further information please read the RAMEN paper (in
preparation).

[`RAMEN::nullDistGE()`](https://ericknavarrod.github.io/RAMEN/reference/nullDistGE.md)
simulates the delta R squared distribution under the null hypothesis of
G and E having no association with DNA methylation (DNAme) variability
through a permutation analysis. To do so, this function shuffles the G
and E variables in the dataset, which is followed by a the variable
selection and modelling steps with selectVariables() and lmGE().These
steps are repeated several times as indicated in the *permutations*
parameter. In other words, by using shuffled G and E data, we simulate
the increase of R2 that would be observed in random data using the RAMEN
methodology.

``` r
# Compute the null distribution 
null_dist = RAMEN::nullDistGE(
  VMRs_df = VMRs_df,
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix = RAMEN::test_environmental_matrix,
  summarized_methyl_VMR = summarized_methyl_VMR,
  permutations = 5,
  covariates = RAMEN::test_covariates,
  seed = 1,
  model_selection = "AIC"
)

#Take a look at the object 
head(null_dist)
#>   VMR_index tot_r_squared model_group R2_difference AIC_difference permutation
#> 1         1     0.2606017           E     0.2495286       157.7547           1
#> 2         2     0.2976822           E     0.2671059       157.5905           1
#> 3         3     0.6086939         G+E     0.5935297       143.5972           1
#> 4         4     0.3777878           E     0.2874781       151.5921           1
#> 5         6     0.6133220         G+E     0.5269145       137.8163           1
#> 6         7     0.3451632           E     0.2497069       155.1677           1
```

The output is a data frame where the most useful column for our purpose
is *R2_difference*, which corresponds to the increase in R squared
obtained by including the G/E variable(s) from the winning model (i.e.,
the R squared difference between the winning model and the model only
with the concomitant variables specified in covariates; tot_r_squared -
basal_rsquared in the lmGE output)

We recommend to use two different thresholds for the winning models
depending of whether they are marginal (G or E) or joint models (G+E or
GxE). The reason for this is that they have different R2_difference
distributions. E and G models have a lower mean R2_difference because
they have a single shuffled term in the model (SNP or E). In comparison,
joint models have a higher mean R2_difference because they have two or
three shuffled terms (SNP, E and SNP\*E), which just by chance increases
their probability of having a higher R2_difference.

``` r
# See the distribution of R2_difference across different winning models
null_dist %>% 
  ggplot2::ggplot(aes(x = R2_difference)) +
  ggplot2::geom_histogram() + 
  ggplot2::facet_grid("model_group") +
  ggplot2::theme_classic()
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![R2 difference (winner - basal) in a suffled data
set.](RAMEN_files/figure-html/unnamed-chunk-13-1.png)

R2 difference (winner - basal) in a suffled data set.

We suggest using the 95th percentile of those distributions as a
threshold to remove bad performing winning models found in our observed
data.

``` r
# Get a cutoff of the 95th percentile of the null distribution for single and joint models
cutoff_single = quantile(null_dist %>% 
                           filter(model_group %in% c("G","E")) %>% 
                           pull(R2_difference), 
                         0.95)
cutoff_joint = quantile(null_dist %>% 
                           filter(model_group %in% c("G+E","GxE")) %>% 
                           pull(R2_difference), 
                         0.95)

#Filter out bad performing winning models.
filtered_res = lmge_res %>% 
  dplyr::mutate(r2_difference_basal = tot_r_squared - basal_rsquared,
         pass_cutoff_threshold = case_when(model_group %in% c("G", "E") ~ r2_difference_basal > cutoff_single,
                                           model_group %in% c("G+E", "GxE") ~ r2_difference_basal > cutoff_joint)) %>% 
  dplyr::filter(pass_cutoff_threshold) %>% #Filter based on the cutoff threshold
  dplyr::select(-pass_cutoff_threshold) #Drop temporary column

#Check the final results
dplyr::glimpse(filtered_res)
#> Rows: 7
#> Columns: 14
#> $ VMR_index           <chr> "5", "6", "28", "39", "96", "98", "106"
#> $ model_group         <chr> "G+E", "GxE", "E", "GxE", "G", "G", "GxE"
#> $ variables           <list> <"21:15671097:T:A", "E43">, <"21:15671097:T:A", "E…
#> $ tot_r_squared       <dbl> 0.7680461, 0.7769276, 0.5555035, 0.7097936, 0.6431…
#> $ g_r_squared         <dbl> 0.39818402, 0.38220942, NA, 0.09101486, 0.6398122…
#> $ e_r_squared         <dbl> 0.2769884, 0.1650993, 0.5244324, 0.2406247, NA, NA…
#> $ gxe_r_squared       <dbl> NA, 0.1432113, NA, 0.3522698, NA, NA, 0.1515714
#> $ AIC                 <dbl> 129.9252, 123.3999, 145.1528, 147.3162, 140.8549, …
#> $ second_winner       <chr> "GxE", "G+E", NA, "G+E", "G+E", NA, "G+E"
#> $ delta_aic           <dbl> 0.07612103, 3.90964213, NA, 3.27939276, 0.34990027…
#> $ delta_r_squared     <dbl> -0.01440811, 0.04856945, NA, 0.05583976, -0.019097…
#> $ basal_AIC           <dbl> 166.8374, 159.6965, 166.5303, 177.6443, 169.6680, …
#> $ basal_rsquared      <dbl> 0.092873702, 0.086407501, 0.031071178, 0.025884145…
#> $ r2_difference_basal <dbl> 0.6751724, 0.6905201, 0.5244324, 0.6839094, 0.6398…
```

We can see that the final data set in this example dropped almost all of
the VMRs we had. This is something expected since we are working with
toy data coming from random sampling.

We recommend the users of the package to include the number of VMRs
where we could not find a conclusive best model in the final results
report (either because no variables were selected with
[`RAMEN::selectVariables()`](https://ericknavarrod.github.io/RAMEN/reference/selectVariables.md)
or because they did not pass the R2_difference threshold obtained with
[`RAMEN::nullDistGE()`](https://ericknavarrod.github.io/RAMEN/reference/nullDistGE.md)).

``` r
# Include the VMRs with no winning model in the final results object
final_res = VMRs_df %>% 
  dplyr::left_join(filtered_res,
                    by = "VMR_index")

# Plot final results
final_res %>% 
  dplyr::group_by(model_group) %>% 
  dplyr::summarise(count = n()) %>% 
  ggplot2::ggplot(aes(x = model_group, y = count)) +
  ggplot2::geom_col() +
  ggplot2::xlab("Best explanatory model") +
  ggplot2::ylab("VMRs") +
  ggplot2::theme_classic() 
```

![Variable Methylated Regions best explanatory
models](RAMEN_files/figure-html/finalresults-1.png)

Variable Methylated Regions best explanatory models

So, we can see that for this toy example, we got the following results:

- VMRs better explained by a G model: 2  
- VMRs better explained by a E model: 1
- VMRs better explained by a G+E model: 1
- VMRs better explained by a GxE model: 3
- VMRs with no conclusive explanatory model: 124

#### Author’s note about model interpretation

For model simplicity, each winning model have a single E and/or G
variable (and its interaction term when applicable). That means that in
a scenario where a given VMR is under the influence of 2 or more Es or
Gs, only the one that better explains the VMR’s DNAme alone will be
selected. In other words, if a VMR in reality is influenced by
e.g. folate intake and smoking, and we have information about both
environmental exposures, the best model (E) will have only folate intake
(in case that is the variable that better explains DNAme variability in
that region alone). So, interpreting this as the VMR not being
potentially under the influence of smoking might not be correct. We
recommend the user to check the total R2 of the winning model to explore
the remaining variance that is not explained by the winning model.

We also stress that **interpretation of individual variables should be
done with caution and used only as exploration and research hypothesis
generation**. Please see [Author’s note about variables
interpretation](#authors-note-about-variables-interpretation)) where we
advice against over-interpretation of the selected variables; the same
logic applies to the variables present in the winning models.

## Variations to the standard workflow

Besides using RAMEN for completing the analysis mentioned above, the
outputs of the package’s function can help users in other DNAme
analyses, such as:

- Reduction of tests prior to an EWAS or differential methylation
  analysis with
  [`RAMEN::findVMRs()`](https://ericknavarrod.github.io/RAMEN/reference/findVMRs.md)(i.e.,
  conducting the analysis on identified VMRs which 1) reduces redundant
  tests by grouping nearby correlated CpGs, and 2) avoids tests in
  non-variant regions). This can help to reduce the multiple hypothesis
  testing burden.
- Summarize a DNAme region of interest with
  [`RAMEN::summarizeVMRs()`](https://ericknavarrod.github.io/RAMEN/reference/summarizeVMRs.md)
- Easily conduct variable selection in high-dimensional data sets to
  identify potentially relevant variables from one or two independent
  data sets with
  [`RAMEN::selectVariables()`](https://ericknavarrod.github.io/RAMEN/reference/selectVariables.md).
- Fit additive and interaction models given two sets of variables of
  interest (not limited to G and E) and select the best explanatory
  model for DNAme data with
  [`RAMEN::selectVariables()`](https://ericknavarrod.github.io/RAMEN/reference/selectVariables.md)
  and
  [`RAMEN::lmGE()`](https://ericknavarrod.github.io/RAMEN/reference/lmGE.md)
  (e.g. exploring the interaction between two environmental dimensions
  and their contribution to DNAme variability).
- Quickly identify SNPs in *cis* of CpG probes with
  [`RAMEN::findCisSNPs()`](https://ericknavarrod.github.io/RAMEN/reference/findCisSNPs.md)
  (e.g. for cis mQTL analyses)
- Get the median correlation of probes in regions of interest (with
  [`RAMEN::medCorVMR()`](https://ericknavarrod.github.io/RAMEN/reference/medCorVMR.md)).

## Frequently Asked Questions

**How can I save the RAMEN data frames that have columns with lists as
observations?**

Saving the data frames produced by RAMEN might seem difficult because it
has lists as observations in several columns, which is not supported by
some read/write functions. We suggest two options:

1.  Use `fwrite()` and `fread()` from the
    *[data.table](https://CRAN.R-project.org/package=data.table)*
    package (recommended because of time and storage space).

``` r
# Example for saving and reading selected_variables object
data.table::fwrite(selected_variables, file = "path/selected_variables.csv")

# Read the csv file and make lists the elements in the required columns
selected_variables = fread("path/selected_variables.csv", data.table = FALSE) %>% 
    mutate(selected_genot = str_split(selected_genot, pattern = "\\|"), # fwrite saves lists as strings separated by |
         selected_env =str_split(selected_env, pattern = "\\|"),
         VMR_index = as.character(VMR_index)) 
```

2.  Save files as .rds

``` r
# Example for saving the selected_variables object
saveRDS(selected_variables, file = "path/selected_variables.Rds")

#Load the object
readRDS(file = "path/selected_variables.Rds")
```

## Session info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] doRNG_1.8.6.2    rngtools_1.5.2   foreach_1.5.2    tidyr_1.3.1     
#> [5] ggplot2_4.0.1    dplyr_1.1.4      RAMEN_1.0.0      knitr_1.50      
#> [9] BiocStyle_2.38.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] sass_0.4.10          generics_0.1.4       shape_1.4.6.1       
#>  [4] stringi_1.8.7        lattice_0.22-7       digest_0.6.39       
#>  [7] magrittr_2.0.4       evaluate_1.0.5       grid_4.5.2          
#> [10] RColorBrewer_1.1-3   bookdown_0.46        iterators_1.0.14    
#> [13] fastmap_1.2.0        Matrix_1.7-4         glmnet_4.1-10       
#> [16] jsonlite_2.0.0       DBI_1.2.3            survival_3.8-3      
#> [19] BiocManager_1.30.27  purrr_1.2.0          scales_1.4.0        
#> [22] relaimpo_2.2-7       codetools_0.2-20     textshaping_1.0.4   
#> [25] jquerylib_0.1.4      cli_3.6.5            mitools_2.4         
#> [28] rlang_1.1.6          splines_4.5.2        withr_3.0.2         
#> [31] cachem_1.1.0         yaml_2.3.12          tools_4.5.2         
#> [34] parallel_4.5.2       corpcor_1.6.10       boot_1.3-32         
#> [37] BiocGenerics_0.56.0  vctrs_0.6.5          R6_2.6.1            
#> [40] stats4_4.5.2         lifecycle_1.0.4      stringr_1.6.0       
#> [43] Seqinfo_1.0.0        S4Vectors_0.48.0     fs_1.6.6            
#> [46] IRanges_2.44.0       MASS_7.3-65          ragg_1.5.0          
#> [49] pkgconfig_2.0.3      desc_1.4.3           pkgdown_2.2.0       
#> [52] pillar_1.11.1        bslib_0.9.0          gtable_0.3.6        
#> [55] Rcpp_1.1.0           glue_1.8.0           systemfonts_1.3.1   
#> [58] xfun_0.54            tibble_3.3.0         GenomicRanges_1.62.1
#> [61] tidyselect_1.2.1     farver_2.1.2         survey_4.4-8        
#> [64] htmltools_0.5.9      rmarkdown_2.30       labeling_0.4.3      
#> [67] compiler_4.5.2       S7_0.2.1
```
