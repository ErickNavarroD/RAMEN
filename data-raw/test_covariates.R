## code to prepare `test_covariates` dataset
set.seed(123)
sample_size = 30
test_covariates = matrix(rnorm(sample_size, 0, 1),
                         nrow = sample_size,
                         ncol = 1)
rownames(test_covariates) = paste("ID", as.character(1:sample_size), sep = "")
colnames(test_covariates) = "covar1"

usethis::use_data(test_covariates, overwrite = TRUE)
