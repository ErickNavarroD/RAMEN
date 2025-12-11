## code to prepare `test_environmental_matrix`
# Simulate environmental data for 100 variables
set.seed(123)
sample_size <- 30
test_environmental_matrix <- matrix(rnorm(100 * sample_size, 0, 1),
  nrow = sample_size,
  ncol = 100
)
rownames(test_environmental_matrix) <- paste("ID",
                                             as.character(1:sample_size),
                                             sep = "")
colnames(test_environmental_matrix) <- paste("E",
                                             as.character(1:100),
                                             sep = "")

usethis::use_data(test_environmental_matrix, overwrite = TRUE)
