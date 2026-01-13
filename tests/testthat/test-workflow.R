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
      var_method = "ultrastable",
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
  expect_error(
    RAMEN::findVML(
      methylation_data = as.matrix(RAMEN::test_methylation_data),
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 0,
      var_method = "variance",
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ),
    "The methylation_data object must be a data frame with samples as columns and probes as rows."
  )
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "a",
      cor_threshold = 0,
      var_method = "variance",
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    )
  )
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 0,
      var_method = "a",
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ),
    "The method must be either 'mad' or 'variance'. Please select one of those options"
  )

})

test_that("findVML works with EPICv2 probes", {
  epic2_methylation_data <- RAMEN::test_methylation_data
  rownames(epic2_methylation_data) <- data.frame(IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations) |>
    dplyr::filter(chr == "chr21") |>
    arrange(chr, pos) |> #Make sure to extract neighbouring probes to have VML
    slice_head(n = nrow(RAMEN::test_methylation_data)) |>
    rownames()

  VML_epic2 <- RAMEN::findVML(
    methylation_data = epic2_methylation_data,
    array_manifest = "IlluminaHumanMethylationEPICv2",
    cor_threshold = 0,
    var_method = "variance",
    var_distribution = "ultrastable",
    var_threshold_percentile = 0.99,
    max_distance = 1000
  )
  expect_true(is.list(VML_epic2))
  expect_true(is.data.frame(VML_epic2$VML))
  expect_equal(ncol(VML_epic2$VML), 10 )
})

test_that("findVML works with var_distribution = 'all' and mad score", {
  VML_allvar <- RAMEN::findVML(
    methylation_data = RAMEN::test_methylation_data,
    array_manifest = "IlluminaHumanMethylationEPICv1",
    cor_threshold = 0,
    var_method = "mad",
    var_distribution = "all",
    var_threshold_percentile = 0.9,
    max_distance = 1000
  )
  expect_true(is.list(VML_allvar))
  expect_true(is.data.frame(VML_allvar$VML))
  expect_equal(ncol(VML_allvar$VML), 10 )
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

test_that("summarizeVML adds VML_index when not present", {
  VML_no_index <- VML$VML |>
    dplyr::select(-VML_index)
  summarized_no_index <- RAMEN::summarizeVML(
    VML_df = VML_no_index,
    methylation_data = test_methylation_data
  )
  expect_true(is.data.frame(summarized_no_index))
  expect_true(all(
    colnames(summarized_no_index) %in% paste0("VML", seq_len(nrow(VML_no_index))))
  )
  expect_equal(nrow(summarized_no_index), ncol(test_methylation_data))
}
)

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
  expect_error(
    RAMEN::summarizeVML(
      VML_df = VML$VML,
      methylation_data = "a"
    ),
    "Please make sure the methylation data is a data frame or matrix with samples as columns and probes as rows."
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

test_that("findCisSNPs adds a VML index when it is not present", {
  VML_cis_snps_noID <- RAMEN::findCisSNPs(
    VML_df = VML$VML |>
      dplyr::select(-VML_index),
    genotype_information = RAMEN::test_genotype_information,
    distance = 1e+06)
  expect_true("VML_index" %in% colnames(VML_cis_snps_noID))
})

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
    RAMEN::findCisSNPs(VML_df = VML$VML,
                       genotype_information = RAMEN::test_genotype_information |>
                         dplyr::select(-CHROM),
                       distance = 1e+06
    ),
    "Please make sure the genotype_information object has the required columns with the appropiate names (check documentation for further information)",
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

test_that("findCisSNPs returns the right number of cis SNPs", {
  VML_test <- data.frame(
    VML_index = "1",
    seqnames = "chr1",
    start = 1000,
    end = 2000,
    type = "VMR"
  )
  genot_info_test <- data.frame(
    CHROM = c("chr1", "chr1", "chr1", "chr1"),
    POS = c(1, 500, 2500, 4000),
    ID = c("rs1", "rs2", "rs3", "rs4")
  )
  test_1 <- RAMEN::findCisSNPs(
    VML_df = VML_test,
    genotype_information = genot_info_test,
    distance = 1
  )
  test_500 <- RAMEN::findCisSNPs(
    VML_df = VML_test,
    genotype_information = genot_info_test,
    distance = 500
  )
  test_1000 <- RAMEN::findCisSNPs(
    VML_df = VML_test,
    genotype_information = genot_info_test,
    distance = 1000
  )
  test_2000 <- RAMEN::findCisSNPs(
    VML_df = VML_test,
    genotype_information = genot_info_test,
    distance = 2000
  )
  expect_equal(test_1$surrounding_SNPs, 0) # no SNPs are within 1 bp
  expect_equal(test_500$surrounding_SNPs, 2) # only rs2 and rs3 are within 500bp
  expect_equal(test_1000$surrounding_SNPs, 3) # rs1, rs2 and rs3 are within 1000bp
  expect_equal(test_2000$surrounding_SNPs, 4) # all 4 snps are within 2000bp
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
test_that("selectVariables output structure is correct", {
  expect_true(is.data.frame(selected_variables))
  expect_equal(ncol(selected_variables), 3)
  expect_equal(nrow(selected_variables), nrow(VML_cis_snps))
  expect_true(all(
    c("VML_index", "selected_genot", "selected_env") %in%
      colnames(selected_variables)
  ))
})

#Test that errors happen when expected
test_that("selectVariables throws errors when expected", {
  expect_error(
    RAMEN::selectVariables(
      VML_df = "a",
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML
    ),
    "Please make sure the VML data frame (VML_df) contains the columns 'SNP' and 'VML_index'.",
    fixed = TRUE
    )
  #Test error when there are ID mismatches
  test_genot <- RAMEN::test_genotype_matrix
  colnames(test_genot) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps,
      genotype_matrix = test_genot,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML
    ),
    "Individual IDs in summarized_methyl_VML do not match individual IDs in genotype_matrix",
    fixed = TRUE
  )
  #Test error when there is argument mismatch with the environmental_matrix
  test_env <- RAMEN::test_environmental_matrix
  rownames(test_env) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = test_env,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML
    ),
    "Individual IDs in summarized_methyl_VML do not match individual IDs in environmental_matrix",
    fixed = TRUE
  )
  #Test error when there is argument mismatch with the covariates
  test_cov <- RAMEN::test_covariates
  rownames(test_cov) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = test_cov,
      summarized_methyl_VML = summarized_methyl_VML
    ),
    "Individual IDs in summarized_methyl_VML do not match individual IDs in the covariates matrix",
    fixed = TRUE
  )

  #Test that matrix arguments throw errors if input is not a matrix
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = as.data.frame(RAMEN::test_environmental_matrix),
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML
    ),
    "Please make sure the environmental data is provided as a matrix.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps,
      genotype_matrix = as.data.frame(RAMEN::test_genotype_matrix),
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML
    ),
    "Please make sure the genotype data is provided as a matrix.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = as.data.frame(RAMEN::test_covariates),
      summarized_methyl_VML = summarized_methyl_VML
    ),
    "Please make sure the covariates data is provided as a matrix.",
    fixed = TRUE
  )
  #Test missing columns in VML_df
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps |>
        dplyr::select(-SNP),
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML
    ),
    "Please make sure the VML data frame (VML_df) contains the columns 'SNP' and 'VML_index'.",
    fixed = TRUE
  )
  #Test missing values in genotype matrix
  #Introduce NA values
  test_genot_na <- RAMEN::test_genotype_matrix
  test_genot_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps,
      genotype_matrix = test_genot_na,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML
    ),
    "Please make sure the genotype matrix contains only finite numeric values.",
    fixed = TRUE
  )
  #Test missing values in environmental matrix
  #Introduce NA values
  test_env_na <- RAMEN::test_environmental_matrix
  test_env_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = test_env_na,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML
    ),
    "Please make sure the environmental matrix contains only finite numeric values.",
    fixed = TRUE
  )
  #Test missing values in covariates matrix
  #Introduce NA values
  test_cov_na <- RAMEN::test_covariates
  test_cov_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = test_cov_na,
      summarized_methyl_VML = summarized_methyl_VML
    ),
    "Please make sure the covariates matrix contains only finite numeric values.",
    fixed = TRUE
  )
  #Test missing values in summarized methylation VML
  test_summeth_na <- summarized_methyl_VML
  test_summeth_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_df = VML_cis_snps,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = test_summeth_na
    ),
    "Please make sure the summarized_methyl_VML data frame contains only finite numeric values.",
    fixed = TRUE
  )
})


#### lmGE() ####
# Use only 10 VML
lmge_res <- RAMEN::lmGE(
  selected_variables = selected_variables[1:10,],
  summarized_methyl_VML = summarized_methyl_VML,
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix = RAMEN::test_environmental_matrix,
  covariates = RAMEN::test_covariates,
  model_selection = "AIC"
)

test_that("lmGE works with BIC model selection", {
  lmge_res_bic <- RAMEN::lmGE(
    selected_variables = selected_variables[1:10,],
    summarized_methyl_VML = summarized_methyl_VML,
    genotype_matrix = RAMEN::test_genotype_matrix,
    environmental_matrix = RAMEN::test_environmental_matrix,
    covariates = RAMEN::test_covariates,
    model_selection = "BIC"
  )
  expect_true(is.data.frame(lmge_res_bic))
  expect_equal(ncol(lmge_res_bic), 13)
  expect_equal(nrow(lmge_res_bic), 10)
})

test_that("lmGE output structure is correct", {
  expect_true(is.data.frame(lmge_res))
  expect_equal(ncol(lmge_res), 13)
  expect_equal(nrow(lmge_res), 10)
})

test_that("lmGE throws errors when expected", {
  expect_error(
    RAMEN::lmGE(
      selected_variables = "a",
      summarized_methyl_VML = summarized_methyl_VML,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      model_selection = "AIC"
    ),
    "Please make sure the selected_variables data frame contains the columns 'VML_index', 'selected_genot' and 'selected_env'.",
    fixed = TRUE
  )
})

#### nullDistGE() ####
permutations <- 2
null_dist <- RAMEN::nullDistGE(
  VML_df = VML_cis_snps[1:10,],
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix = RAMEN::test_environmental_matrix,
  summarized_methyl_VML = summarized_methyl_VML,
  permutations = permutations,
  covariates = RAMEN::test_covariates,
  seed = 1,
  model_selection = "AIC"
)
test_that("nullDistGE output structure is correct", {
  expect_true(is.data.frame(null_dist))
  expect_equal(ncol(null_dist), 6)
  expect_equal(nrow(null_dist), nrow(VML_cis_snps[1:10,])*permutations)
})

test_that("nullDistGE works with BIC", {
  null_dist_bic <- RAMEN::nullDistGE(
    VML_df = VML_cis_snps[1:10,],
    genotype_matrix = RAMEN::test_genotype_matrix,
    environmental_matrix = RAMEN::test_environmental_matrix,
    summarized_methyl_VML = summarized_methyl_VML,
    permutations = permutations,
    covariates = RAMEN::test_covariates,
    seed = 1,
    model_selection = "BIC"
  )
  expect_true(is.data.frame(null_dist_bic))
  expect_equal(ncol(null_dist_bic), 6)
  expect_equal(nrow(null_dist_bic), nrow(VML_cis_snps[1:10,])*permutations)
})

test_that("nullDistGE throws errors when expected", {
  expect_error(
    RAMEN::nullDistGE(
      VML_df = "a",
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      summarized_methyl_VML = summarized_methyl_VML,
      permutations = 2,
      covariates = RAMEN::test_covariates,
      seed = 1,
      model_selection = "AIC"
    ),
    "Please make sure the VML data frame (VML_df) contains the columns 'SNP' and 'VML_index'.",
    fixed = TRUE
  )
  #Test error when genotype_matrix has NA
  #Introduce NA values
  test_genot_na <- RAMEN::test_genotype_matrix
  test_genot_na[1, 1] <- NA
  expect_error(
    RAMEN::nullDistGE(
      VML_df = VML_cis_snps[1:10,],
      genotype_matrix = test_genot_na,
      environmental_matrix = RAMEN::test_environmental_matrix,
      summarized_methyl_VML = summarized_methyl_VML,
      permutations = 2,
      covariates = RAMEN::test_covariates,
      seed = 1,
      model_selection = "AIC"
    ),
    "Please make sure the genotype matrix contains only finite numeric values.",
    fixed = TRUE
  )
  #Test error when environmental_matrix has NA
  #Introduce NA values
  test_env_na <- RAMEN::test_environmental_matrix
  test_env_na[1, 1] <- NA
  expect_error(
    RAMEN::nullDistGE(
      VML_df = VML_cis_snps[1:10,],
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = test_env_na,
      summarized_methyl_VML = summarized_methyl_VML,
      permutations = 2,
      covariates = RAMEN::test_covariates,
      seed = 1,
      model_selection = "AIC"
    ),
    "Please make sure the environmental matrix contains only finite numeric values.",
    fixed = TRUE
  )
  #Test error when covariates has NA
  #Introduce NA values
  test_cov_na <- RAMEN::test_covariates
  test_cov_na[1, 1] <- NA
  expect_error(
    RAMEN::nullDistGE(
      VML_df = VML_cis_snps[1:10,],
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      summarized_methyl_VML = summarized_methyl_VML,
      permutations = 2,
      covariates = test_cov_na,
      seed = 1,
      model_selection = "AIC"
    ),
    "Please make sure the covariates matrix contains only finite numeric values.",
    fixed = TRUE
  )
  #Test error when summarized_methyl_VML has NA
  #Introduce NA values
  test_summeth_na <- summarized_methyl_VML
  test_summeth_na[1, 1] <- NA
  expect_error(
    RAMEN::nullDistGE(
      VML_df = VML_cis_snps[1:10,],
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      summarized_methyl_VML = test_summeth_na,
      permutations = 2,
      covariates = RAMEN::test_covariates,
      seed = 1,
      model_selection = "AIC"
    ),
    "Please make sure the summarized_methyl_VML data frame contains only finite numeric values.",
    fixed = TRUE
  )
}
)

#### Clean environment ####
rm(VML, lmge_res, null_dist, summarized_methyl_VML, selected_variables, VML_cis_snps, permutations)
