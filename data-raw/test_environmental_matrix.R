## code to prepare `test_environmental_matrix`
# Simulate environmental data for 100 variables
set.seed(123)
sample_size <- 30
test_environmental_matrix <- matrix(rnorm(100 * sample_size, 0, 1),
  nrow = sample_size,
  ncol = 100
)
rownames(test_environmental_matrix) <- paste0("ID",
  as.character(1:sample_size)
)
colnames(test_environmental_matrix) <- paste0("E",
  as.character(1:100)
)

usethis::use_data(test_environmental_matrix, overwrite = TRUE)
