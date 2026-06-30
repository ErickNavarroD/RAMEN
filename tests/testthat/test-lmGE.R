test_that("lmGE output structure is correct", {
  lmge_res <- RAMEN::lmGE(
    selected_variables = selected_variables_test[1:10, ],
    summarized_methyl_VML = summarized_methyl_VML_test,
    genotype_matrix = RAMEN::test_genotype_matrix,
    environmental_matrix = RAMEN::test_environmental_matrix,
    covariates = RAMEN::test_covariates,
    model_selection = "AIC"
  )
  expect_true(is.data.frame(lmge_res))
  expect_equal(ncol(lmge_res), 13)
  expect_equal(nrow(lmge_res), 10)
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
    "Please make sure the selected_variables data frame contains the columns 'VML_index', 'selected_genot' and 'selected_env'.",
    fixed = TRUE
  )
})
