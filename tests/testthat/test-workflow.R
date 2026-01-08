# Since the input of some functions is the output of others in the package,
# the whole workflow is going to be tested in this script to minimize the
# run time.

library(testthat)
library(dplyr)

#### findVML() ####
VML <- RAMEN::findVML(
  methylation_data = RAMEN::test_methylation_data,
  array_manifest = "IlluminaHumanMethylationEPICv1",
  cor_threshold = 0,
  var_method = "variance",
  var_distribution = "ultrastable",
  var_threshold_percentile = 0.99,
  max_distance = 1000
)

test_that("findVML variance calculation is correct", {
  probe_test <- VML$highly_variable_probes[1:10, ]
  observed_variance <- RAMEN::test_methylation_data[probe_test$TargetID, ] |>
    apply(1, var)
  names(observed_variance) <- NULL
  expect_equal(observed_variance, probe_test$var_score)
})

test_that("findVML output structure is correct", {
  expect_true(is.list(VML))
  expect_true("VML" %in% names(VML))
  expect_true("highly_variable_probes" %in% names(VML))
  expect_true("var_score_threshold" %in% names(VML))
  expect_true(is.data.frame(VML$VML))
  expect_true(is.data.frame(VML$highly_variable_probes))
})

test_that("findVML handles a different var_method option", {
  VML_result_mad <- RAMEN::findVML(
    methylation_data = RAMEN::test_methylation_data,
    array_manifest = "IlluminaHumanMethylationEPICv1",
    cor_threshold = 0,
    var_method = "mad",
    var_distribution = "ultrastable",
    var_threshold_percentile = 0.99,
    max_distance = 1000
  )
  probe_test <- VML_result_mad$highly_variable_probes[1:10, ]
  observed_mad <- RAMEN::test_methylation_data[probe_test$TargetID, ] |>
    apply(1, mad)
  names(observed_mad) <- NULL
  expect_equal(observed_mad, probe_test$var_score)
  expect_true(is.list(VML_result_mad))
  expect_true("VML" %in% names(VML_result_mad))
})

test_that("findVML throws errors when expected", {
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 1000,
      var_method = "invalid_method",
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ),
    "'cor_threshold' must be of type 'numeric' and from 0 to 1"
  )
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 0,
      var_method = "variance",
      var_distribution = "binomial",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ),
    "'var_distribution' must be one of 'all' or 'ultrastable'"
  )
})

test_that("sVMPs have no correlation", {
  sVMPs <- VML$VML |>
    dplyr::filter(type == "sVMP") |>
    dplyr::pull(median_correlation)
  expect_true(all(is.na(sVMPs)))
})

test_that("correlation is computed correctly", {
  VML_test <- VML$VML |>
    dplyr::filter(type == "VMR") |>
    dplyr::arrange(n_VMPs) |>
    dplyr::slice_tail(n = 1) # Get the VMR with highest number of HVPs
  probes <- unlist(VML_test$probes)
  methylation_subset <- RAMEN::test_methylation_data[probes, ]
  cor_matrix <- cor(t(methylation_subset))
  cor_values <- cor_matrix[lower.tri(cor_matrix)]
  median_cor <- median(cor_values)
  expect_equal(VML_test$median_correlation, median_cor)
})

#### summarizeVML() ####
summarized_methyl_VML <- RAMEN::summarizeVML(
  VML_df = VML$VML,
  methylation_data = test_methylation_data
)

test_that("summarizeVML output structure is correct", {
  expect_true(is.data.frame(summarized_methyl_VML))
  expect_equal(ncol(summarized_methyl_VML), nrow(VML$VML))
  expect_equal(nrow(summarized_methyl_VML), ncol(test_methylation_data))
})

test_that("summarizeVML values are correct", {
  # First for sVMPs: the summarized value should be equal to the methylation
  VML_test <- VML$VML |>
    dplyr::filter(type == "sVMP") |>
    dplyr::slice_head(n = 1) # Get the first sVMP
  probe <- unlist(VML_test$probes)
  expected <- RAMEN::test_methylation_data[probe, ] |> unlist()
  observed <- summarized_methyl_VML[, VML_test$VML_index]
  names(observed) <- rownames(summarized_methyl_VML)
  expect_equal(observed, expected)
  # now for VMRs: the summarized value should be the median across probes
  VMR_test <- VML$VML |>
    dplyr::filter(type == "VMR") |>
    dplyr::slice_head(n = 1) # Get the first VMR
  probes <- unlist(VMR_test$probes)
  expected <- apply(
    RAMEN::test_methylation_data[probes, ],
    2,
    median
  )
  observed <- summarized_methyl_VML[, VMR_test$VML_index]
  names(observed) <- rownames(summarized_methyl_VML)
  expect_equal(observed, expected)
})

test_that("summarizeVML throws errors when expected", {
  expect_error(
    RAMEN::summarizeVML(
      VML_df = "a",
      methylation_data = test_methylation_data
    ),
    "Please provide a data frame in VML_df"
  )
})

test_that("summarizeVML works when methylation_data is a matrix", {
  summarized_methyl_VML_matrix <- RAMEN::summarizeVML(
    VML_df = VML$VML,
    methylation_data = as.matrix(test_methylation_data)
  )
  expect_true(is.data.frame(summarized_methyl_VML_matrix))
  expect_equal(ncol(summarized_methyl_VML_matrix), nrow(VML$VML))
  expect_equal(nrow(summarized_methyl_VML_matrix), ncol(test_methylation_data))
  svmp <- VML$VML |>
    dplyr::filter(type == "sVMP") |>
    dplyr::slice_head(n = 1) # Get the first sVMP
  expect_equal(
    summarized_methyl_VML_matrix[, svmp$VML_index],
    test_methylation_data[unlist(svmp$probes), ] |> unlist() |> unname()
  )
})

#### findCisSNPs() ####
VML_cis_snps <- RAMEN::findCisSNPs(
  VML_df = VML$VML,
  genotype_information = RAMEN::test_genotype_information,
  distance = 1e+06
)

test_that("findCisSNPs output structure is correct", {
  expect_true(is.data.frame(VML_cis_snps))
  expect_equal(ncol(VML_cis_snps), ncol(VML$VML) + 2)
  expect_equal(nrow(VML_cis_snps), nrow(VML$VML))
  expect_true(all(
    c(colnames(VML$VML), "surrounding_SNPs", "SNP") %in%
      colnames(VML_cis_snps)
  ))
})

test_that("findCisSNPs throws errors when expected", {
  expect_error(
    RAMEN::findCisSNPs(VML_df = VML$VML |>
                         dplyr::select(-seqnames),
                       genotype_information = RAMEN::test_genotype_information,
                       distance = 1e+06
                       ),
    "Please make sure the VML_df object has the required columns with the appropiate names (check documentation for further information)",
    fixed = TRUE
  )
  expect_error(
    RAMEN::findCisSNPs(VML_df = "a",
                       genotype_information = RAMEN::test_genotype_information,
                       distance = 1e+06
    ),
    "Please make sure the VML_df object is a data frame.",
    fixed = TRUE)
  expect_error(
    RAMEN::findCisSNPs(VML_df = VML$VML,
                       genotype_information = "a",
                       distance = 1e+06
    ),
    "Please make sure the genotype_information object is a data frame.")
}
)

#### selectVariables() ####
selected_variables <- RAMEN::selectVariables(
  VML_df = VML_cis_snps,
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix = RAMEN::test_environmental_matrix,
  covariates = RAMEN::test_covariates,
  summarized_methyl_VML = summarized_methyl_VML,
  seed = 1
)
# test_that("selectVariables output structure is correct", {
#   expect_true(is.data.frame(selected_variables))
#   expect_equal(ncol(selected_variables), 4)
#   expect_equal(nrow(selected_variables), nrow(VML_cis_snps))
#   expect_true(all(
#     c("VML_index", "SNP", "environmental_variable", "model_type") %in%
#       colnames(selected_variables)
#   ))
# })
# Test that a matrix with NA throws an error
# test_that("selectVariables throws errors when expected", {
#   expect_error(
#     RAMEN::selectVariables(
#       VML_df = VML_cis_snps,
#       genotype_matrix = RAMEN::test_genotype_matrix,
#       environmental_matrix = RAMEN::test_environmental_matrix,
#       covariates = RAMEN::test_covariates,
#       summarized_methyl_VML = summarized_methyl_VML |>
#         as.matrix() |>
#         {.[1,1] <- NA; .}, # introduce an NA
#       seed = 1
#     ),
#     "summarized_methyl_VML contains NA values. Please remove or impute them before proceeding."
#   )
# })

#### lmGE() ####
lmge_res <- RAMEN::lmGE(
  selected_variables = selected_variables,
  summarized_methyl_VML = summarized_methyl_VML,
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix = RAMEN::test_environmental_matrix,
  covariates = RAMEN::test_covariates,
  model_selection = "AIC"
)

#### nullDistGE() ####
null_dist <- RAMEN::nullDistGE(
  VML_df = VML_cis_snps,
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix = RAMEN::test_environmental_matrix,
  summarized_methyl_VML = summarized_methyl_VML,
  permutations = 2,
  covariates = RAMEN::test_covariates,
  seed = 1,
  model_selection = "AIC"
)

#### Clean environment ####
