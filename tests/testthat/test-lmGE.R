test_that("lmGE output structure is correct", {
  # Set the parallel backend to use 2 workers
  doParallel::registerDoParallel(2)
  lmge_res <- RAMEN::lmGE(
    selected_variables = selected_variables_test[1:5, ],
    summarized_methyl_VML = summarized_methyl_VML_test,
    genotype_matrix = RAMEN::test_genotype_matrix,
    environmental_matrix = RAMEN::test_environmental_matrix,
    covariates = RAMEN::test_covariates,
    model_selection = "AIC"
  )
  expect_true(is.data.frame(lmge_res))
  expect_equal(ncol(lmge_res), 13)
  expect_equal(nrow(lmge_res), 5)
})

test_that("lmGE throws errors when expected", {
  expect_error(
    RAMEN::lmGE(
      selected_variables = "a",
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      model_selection = "AIC"
    ),
    "Please make sure the input selected_variables belongs to the data.frame class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = "summarized_methyl_VML_test",
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      model_selection = "AIC"
    ),
    "Please make sure the input summarized_methyl_VML belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = "test_genotype_matrix",
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      model_selection = "AIC"
    ),
    "Please make sure the input genotype_matrix belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = "test_environmental_matrix",
      covariates = RAMEN::test_covariates,
      model_selection = "AIC"
    ),
    "Please make sure the input environmental_matrix belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = "test_covariates",
      model_selection = "AIC"
    ),
    "Please make sure the input covariates belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      model_selection = "a"
    ),
    "Please make sure the input model_selection is one of the following options: AIC, BIC .",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      model_selection = 1
    ),
    "Please make sure the input model_selection belongs to the character class.",
    fixed = TRUE
  )
})
