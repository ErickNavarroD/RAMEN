test_that("selectVariables output structure is correct", {
  expect_true(is.data.frame(selected_variables_test))
  expect_equal(ncol(selected_variables_test), 3)
  expect_equal(nrow(selected_variables_test), length(VML_cis_snps_test))
  expect_true(all(
    c("VML_index", "selected_genot", "selected_env") %in%
      colnames(selected_variables_test)
  ))
})

# Test that errors happen when expected
test_that("selectVariables throws errors when expected", {
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = "a",
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the input VML_wSNPs belongs to the GRanges class.",
    fixed = TRUE
  )
  # Test error when there are ID mismatches
  test_genot <- RAMEN::test_genotype_matrix
  colnames(test_genot) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = test_genot,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "The objects rownames(summarized_methyl_VML) and colnames(genotype_matrix) must match.",
    fixed = TRUE
  )
  # Test error when there is argument mismatch with the environmental_matrix
  test_env <- RAMEN::test_environmental_matrix
  rownames(test_env) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = test_env,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "The objects rownames(summarized_methyl_VML) and rownames(environmental_matrix) must match.",
    fixed = TRUE
  )
  # Test error when there is argument mismatch with the covariates
  test_cov <- RAMEN::test_covariates
  rownames(test_cov) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = test_cov,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "The objects rownames(summarized_methyl_VML) and rownames(covariates) must match.",
    fixed = TRUE
  )

  # Test that matrix arguments throw errors if input is not a matrix
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = as.data.frame(RAMEN::test_environmental_matrix),
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the input environmental_matrix belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = as.data.frame(RAMEN::test_genotype_matrix),
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the input genotype_matrix belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = as.data.frame(RAMEN::test_covariates),
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the input covariates belongs to the matrix class.",
    fixed = TRUE
  )
  # Test missing columns in VML_df
  VML_no_SNP = VML_cis_snps_test
  VML_no_SNP$SNP <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_no_SNP,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "The object S4Vectors::mcols(VML_wSNPs) does not have the required columns: VML_index, SNP .",
    fixed = TRUE
  )
  # Test missing values in genotype matrix
  # Introduce NA values
  test_genot_na <- RAMEN::test_genotype_matrix
  test_genot_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = test_genot_na,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the object genotype_matrix contains only finite numeric values (i.e., no NA, NaN or Inf)",
    fixed = TRUE
  )
  # Test missing values in environmental matrix
  # Introduce Inf value
  test_env_na <- RAMEN::test_environmental_matrix
  test_env_na[1, 1] <- Inf
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = test_env_na,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the object environmental_matrix contains only finite numeric values (i.e., no NA, NaN or Inf)",
    fixed = TRUE
  )
  # Test missing values in covariates matrix
  # Introduce NA values
  test_cov_na <- RAMEN::test_covariates
  test_cov_na[1, 1] <- NaN
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = test_cov_na,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the object covariates contains only finite numeric values (i.e., no NA, NaN or Inf)",
    fixed = TRUE
  )
  # Test missing values in summarized methylation VML
  test_summeth_na <- summarized_methyl_VML_test
  test_summeth_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = test_summeth_na
    ),
    "Please make sure the object summarized_methyl_VML contains only finite numeric values (i.e., no NA, NaN or Inf)",
    fixed = TRUE
  )
})
