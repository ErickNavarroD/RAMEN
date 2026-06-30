test_that("nullDistGE output structure is correct", {
  null_dist <- RAMEN::nullDistGE(
    VML_df = VML_cis_snps_test[1:10, ],
    genotype_matrix = RAMEN::test_genotype_matrix,
    environmental_matrix = RAMEN::test_environmental_matrix,
    summarized_methyl_VML = summarized_methyl_VML_test,
    permutations = 2,
    covariates = RAMEN::test_covariates,
    seed = 1,
    model_selection = "AIC"
  ) |>
    suppressMessages()
  expect_true(is.data.frame(null_dist))
  expect_equal(ncol(null_dist), 6)
  expect_equal(nrow(null_dist), nrow(VML_cis_snps_test[1:10, ]) * 2)
})

test_that("nullDistGE works with BIC", {
  null_dist_bic <- RAMEN::nullDistGE(
    VML_df = VML_cis_snps_test[1:10, ],
    genotype_matrix = RAMEN::test_genotype_matrix,
    environmental_matrix = RAMEN::test_environmental_matrix,
    summarized_methyl_VML = summarized_methyl_VML_test,
    permutations = 2,
    covariates = RAMEN::test_covariates,
    seed = 1,
    model_selection = "BIC"
  ) |>
    suppressMessages()
  expect_true(is.data.frame(null_dist_bic))
  expect_equal(ncol(null_dist_bic), 6)
  expect_equal(nrow(null_dist_bic), nrow(VML_cis_snps_test[1:10, ]) * 2)
})

test_that("nullDistGE throws errors when expected", {
  expect_error(
    RAMEN::nullDistGE(
      VML_df = "a",
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      summarized_methyl_VML = summarized_methyl_VML_test,
      permutations = 2,
      covariates = RAMEN::test_covariates,
      seed = 1,
      model_selection = "AIC"
    ) |>
      suppressMessages(),
    "Please make sure the VML data frame (VML_df) contains the columns 'SNP' and 'VML_index'.",
    fixed = TRUE
  )
  # Test error when genotype_matrix has NA
  # Introduce NA values
  test_genot_na <- RAMEN::test_genotype_matrix
  test_genot_na[1, 1] <- NA
  expect_error(
    RAMEN::nullDistGE(
      VML_df = VML_cis_snps_test[1:10, ],
      genotype_matrix = test_genot_na,
      environmental_matrix = RAMEN::test_environmental_matrix,
      summarized_methyl_VML = summarized_methyl_VML_test,
      permutations = 2,
      covariates = RAMEN::test_covariates,
      seed = 1,
      model_selection = "AIC"
    ) |>
      suppressMessages(),
    "Please make sure the genotype matrix contains only finite numeric values.",
    fixed = TRUE
  )
  # Test error when environmental_matrix has NA
  # Introduce NA values
  test_env_na <- RAMEN::test_environmental_matrix
  test_env_na[1, 1] <- NA
  expect_error(
    RAMEN::nullDistGE(
      VML_df = VML_cis_snps_test[1:10, ],
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = test_env_na,
      summarized_methyl_VML = summarized_methyl_VML_test,
      permutations = 2,
      covariates = RAMEN::test_covariates,
      seed = 1,
      model_selection = "AIC"
    ) |>
      suppressMessages(),
    "Please make sure the environmental matrix contains only finite numeric values.",
    fixed = TRUE
  )
  # Test error when covariates has NA
  # Introduce NA values
  test_cov_na <- RAMEN::test_covariates
  test_cov_na[1, 1] <- NA
  expect_error(
    RAMEN::nullDistGE(
      VML_df = VML_cis_snps_test[1:10, ],
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      summarized_methyl_VML = summarized_methyl_VML_test,
      permutations = 2,
      covariates = test_cov_na,
      seed = 1,
      model_selection = "AIC"
    ) |>
      suppressMessages(),
    "Please make sure the covariates matrix contains only finite numeric values.",
    fixed = TRUE
  )
  # Test error when summarized_methyl_VML has NA
  # Introduce NA values
  test_summeth_na <- summarized_methyl_VML_test
  test_summeth_na[1, 1] <- NA
  expect_error(
    RAMEN::nullDistGE(
      VML_df = VML_cis_snps_test[1:10, ],
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      summarized_methyl_VML = test_summeth_na,
      permutations = 2,
      covariates = RAMEN::test_covariates,
      seed = 1,
      model_selection = "AIC"
    ) |>
      suppressMessages(),
    "Please make sure the summarized_methyl_VML data frame contains only finite numeric values.",
    fixed = TRUE
  )
})
