# Changelog

## RAMEN 2.0.0

In this version, we have made an important change in RAMEN terminology
across all the code and documentation to more accurately reflect the
biological concepts represented by the data. The term “Variably
Methylated Regions (VMR)” used in RAMEN v1 has been replaced by
“Variably Methylated Loci (VML)” in RAMEN v2, as not all VML are
composed of 2 or more highly variable probes. VML are further composed
of Variably Methylated Regions (previously named “canonical VMR” in
RAMEN v1) and sparse Variably Methylated Probes (sVMPs; previously named
“non-canonical VMR” in RAMENv1). To be clear, there are no changes in
how these VML are identified, we only changed how we label these
categories.

| Updated name in RAMEN v2                | Deprecated name in RAMEN v1      |
|-----------------------------------------|----------------------------------|
| Variably Methylated Loci (VML)          | Variably Methylated Region (VMR) |
| Variably Methylated Region (VMR)        | canonical VMR (cVMR)             |
| sparse Variably Methylated Probe (sVMP) | non-canonical VMR (ncVMR)        |

Terminology update

- To reflect the terminology change, the following functions had a name
  change:
  [`findVML()`](https://ericknavarrod.github.io/RAMEN/reference/findVML.md)
  (previously named `findVMRs()` in RAMEN v1) and
  [`summarizeVML()`](https://ericknavarrod.github.io/RAMEN/reference/summarizeVML.md)
  (previously named `summarizeVMRs()` in RAMEN v1).

- [`findVML()`](https://ericknavarrod.github.io/RAMEN/reference/findVML.md):

  - Output: list does not separate VMRs and sVMPs into two different
    list elements anymore. Now, a single element (“VML”) is returned in
    the output list, which contains both VMRs and sVMPs, labelled
    accordingly under the *type* column; this VML element is now a data
    frame, and not a Genomic Ranges object to facilitate data wrangling
    and plotting. The function now automatically indexes the VML.

  - The user does not need to provide the array manifest anymore if
    working with the Illumina 450k, EPICv1 or EPICv2 array. The
    `array_manifest` argument accepts now
    “IlluminaHumanMethylation450k”, “IlluminaHumanMethylationEPICv1” and
    “IlluminaHumanMethylationEPICv2”.

  - There is a new method to identify VML using ultrastable probes
    (probes which DNA methylation is known to be stable independently of
    tissue and developmental stage) to discriminate Highly Variable
    Probes, which are then grouped into VML. This method is the default
    one now. For more information please see the
    [`findVML()`](https://ericknavarrod.github.io/RAMEN/reference/findVML.md)
    documentation and the package vignette. The previously default
    method to identify Highly Variable Probes (top 10% of probes with
    the highest variance in the data set) is still available using the
    argument `var_distribution = "all"`.

- [`nullDistGE()`](https://ericknavarrod.github.io/RAMEN/reference/nullDistGE.md):
  Prints messages to keep track of the progress. Fixed a bug that made
  doFuture parallelization strategies crash.

- All functions have examples in the documentation.

- Added tests to reach a code coverage of \>90% in all functions.

- Improved error catches to make functions stop early when the inputs
  are not in the right format. Fixed various bugs throughout the code
  (no user.

- Added news, citation and contributing files to the repository.

- Citation info is provided when loading the package.

- The package repository has now informative badges and Continuous
  Integration checks.
