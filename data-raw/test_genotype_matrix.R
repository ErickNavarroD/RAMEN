## code to prepare `test_genotype_matrix` dataset goes here
# This code makes use of test_genotype_information, which was created by just
# extracting the SNP positions from a real private data set to have an example
# of the SNP IDs.
load(test_genotype_information.Rdata)

set.seed(123)
test_genotype_matrix <- matrix(
  rbinom(nrow(test_genotype_information) * sample_size,
         2,
         0.5
         ),
  ncol = sample_size,
  nrow = nrow(test_genotype_information)
)
colnames(test_genotype_matrix) <- paste("ID",
                                        as.character(1:sample_size),
                                        sep = ""
                                        )
rownames(test_genotype_matrix) <- test_genotype_information$ID

usethis::use_data(test_genotype_matrix, overwrite = TRUE)
