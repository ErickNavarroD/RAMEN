#' Genotype matrix example
#'
#' Genotype matrix example using a gene-dosage model, which encodes the SNPs ordinally depending on the genotype allele charge, such as 2 (AA), 1 (AB) and 0 (BB). Valus were drawn from a binomial distribution with size 2 and probability 0.5.
#'
#' @format ## `test_genotype_matrix`
#' A data frame with 8,539 rows and 30 columns:
#' \describe{
#'   \item{*rownames*}{SNP ID}
#'   \item{ID1:30}{Individual's 1 to 30 genotypes; column names correspond to individual IDs}
#'   ...
#' }
"test_genotype_matrix"
