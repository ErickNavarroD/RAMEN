% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_genotype_matrix.R
\docType{data}
\name{test_genotype_matrix}
\alias{test_genotype_matrix}
\title{Genotype matrix example}
\format{
\subsection{\code{test_genotype_matrix}}{

A data frame with 8,539 rows and 30 columns:
\describe{
\item{\emph{rownames}}{SNP ID}
\item{ID1:30}{Individual's 1 to 30 genotypes; column names correspond to individual IDs}
...
}
}
}
\usage{
test_genotype_matrix
}
\description{
Genotype matrix example using a gene-dosage model, which encodes the SNPs ordinally depending on the genotype allele charge, such as 2 (AA), 1 (AB) and 0 (BB). Valus were drawn from a binomial distribution with size 2 and probability 0.5.
}
\keyword{datasets}
