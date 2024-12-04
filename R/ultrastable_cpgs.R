#' Ultrastable probes
#'
#' This data set contains the list of ultrastable probes identified by [Rachel Edgar et. al.,(2014)](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/1756-8935-7-28). This publication identified ultrastable CpGs across many tissues and conditions using the Illumina 450k array. Ultrastable probes are defined as CpGs consistently methylated or unmethylated in every sample (1,737 samples from 30 publically available studies). These CpGs are used to create a "null DNAme variance" distribution in the RAMEN package, from which a threshold is taken to identify Highly Variable Probes.
#'
#' @format ## `ultrastable_cpgs`
#' A vector with the name of the 15,224 ultrastable probes identified by Edgar et al. (2014). The name of the probes are based on the Illumina 450k manifest.
#'
#' @source https://static-content.springer.com/esm/art%3A10.1186%2F1756-8935-7-28/MediaObjects/13072_2014_333_MOESM2_ESM.txt
"ultrastable_cpgs"
