test_that("selectVariables output structure is correct", {
  expect_true(is.data.frame(selected_variables_test))
  expect_equal(ncol(selected_variables_test), 3)
  expect_equal(nrow(selected_variables_test), nrow(VML_cis_snps_test))
  expect_true(all(
    c("VML_index", "selected_genot", "selected_env") %in%
      colnames(selected_variables_test)
  ))
})

# Test that errors happen when expected
test_that("selectVariables throws errors when expected", {
  expect_error(
    RAMEN::selectVariables(
      VML_df = "a",
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the VML data frame (VML_df) contains the columns 'SNP' and 'VML_index'.",
    fixed = TRUE
  )
  # Test error when there are ID mismatches
  test_genot <- RAMEN::test_genotype_matrix
  colnames(test_genot) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps_test,
      genotype_matrix = test_genot,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Individual IDs in summarized_methyl_VML do not match individual IDs in genotype_matrix",
    fixed = TRUE
  )
  # Test error when there is argument mismatch with the environmental_matrix
  test_env <- RAMEN::test_environmental_matrix
  rownames(test_env) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = test_env,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Individual IDs in summarized_methyl_VML do not match individual IDs in environmental_matrix",
    fixed = TRUE
  )
  # Test error when there is argument mismatch with the covariates
  test_cov <- RAMEN::test_covariates
  rownames(test_cov) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = test_cov,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Individual IDs in summarized_methyl_VML do not match individual IDs in the covariates matrix",
    fixed = TRUE
  )

  # Test that matrix arguments throw errors if input is not a matrix
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = as.data.frame(RAMEN::test_environmental_matrix),
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the environmental data is provided as a matrix.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps_test,
      genotype_matrix = as.data.frame(RAMEN::test_genotype_matrix),
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the genotype data is provided as a matrix.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = as.data.frame(RAMEN::test_covariates),
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the covariates data is provided as a matrix.",
    fixed = TRUE
  )
  # Test missing columns in VML_df
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps_test |>
        dplyr::select(-SNP),
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the VML data frame (VML_df) contains the columns 'SNP' and 'VML_index'.",
    fixed = TRUE
  )
  # Test missing values in genotype matrix
  # Introduce NA values
  test_genot_na <- RAMEN::test_genotype_matrix
  test_genot_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps_test,
      genotype_matrix = test_genot_na,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the genotype matrix contains only finite numeric values.",
    fixed = TRUE
  )
  # Test missing values in environmental matrix
  # Introduce NA values
  test_env_na <- RAMEN::test_environmental_matrix
  test_env_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = test_env_na,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the environmental matrix contains only finite numeric values.",
    fixed = TRUE
  )
  # Test missing values in covariates matrix
  # Introduce NA values
  test_cov_na <- RAMEN::test_covariates
  test_cov_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = test_cov_na,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the covariates matrix contains only finite numeric values.",
    fixed = TRUE
  )
  # Test missing values in summarized methylation VML
  test_summeth_na <- summarized_methyl_VML_test
  test_summeth_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = test_summeth_na
    ),
    "Please make sure the summarized_methyl_VML data frame contains only finite numeric values.",
    fixed = TRUE
  )
})
