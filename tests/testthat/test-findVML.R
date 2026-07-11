test_that("findVML variance calculation is correct", {
  probe_test <- VML_test$highly_variable_probes[1:10, ]
  observed_variance <- RAMEN::test_methylation_data[probe_test$TargetID, ] |>
    apply(1, var)
  names(observed_variance) <- NULL
  expect_equal(observed_variance, probe_test$var_score)
})

test_that("findVML output structure is correct", {
  expect_true(is.list(VML_test))
  expect_true("VML" %in% names(VML_test))
  expect_true("highly_variable_probes" %in% names(VML_test))
  expect_true("var_score_threshold" %in% names(VML_test))
  expect_true(is(VML_test$VML, "GRanges"))
  expect_true(is.data.frame(VML_test$highly_variable_probes))
})

test_that("findVML handles a different var_method option", {
  # Set the parallel backend to use 2 workers
  doParallel::registerDoParallel(2)
  VML_result_mad <- RAMEN::findVML(
    methylation_data = RAMEN::test_methylation_data,
    array_manifest = "IlluminaHumanMethylationEPICv1",
    cor_threshold = 0,
    var_method = "mad",
    var_distribution = "ultrastable",
    var_threshold_percentile = 0.99,
    max_distance = 1000
  ) |>
    suppressMessages()
  probe_test <- VML_result_mad$highly_variable_probes[1:10, ]
  observed_mad <- RAMEN::test_methylation_data[probe_test$TargetID, ] |>
    apply(1, mad)
  names(observed_mad) <- NULL
  expect_equal(observed_mad, probe_test$var_score)
  expect_true(is.list(VML_result_mad))
  expect_true("VML" %in% names(VML_result_mad))
})

test_that("findVML throws errors when expected", {
  #### methylation_data ####
  expect_error(
    RAMEN::findVML(
      methylation_data = as.matrix(RAMEN::test_methylation_data),
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 0,
      var_method = "variance",
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ) |>
      suppressMessages(),
    "Please make sure the input methylation_data belongs to the data.frame class."
  )
  #### array_manifest ####
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "a",
      cor_threshold = 0,
      var_method = "variance",
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ) |>
      suppressMessages(),
    paste("Please make sure the input array_manifest is one of the following",
          "options: IlluminaHumanMethylation450k, IlluminaHumanMethylationEPICv1,",
          "IlluminaHumanMethylationEPICv2 . Otherwise, please provide a custom",
          "manifest as a data frame."),
    fixed = TRUE
  )
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = data.frame(),
      cor_threshold = 0,
      var_method = "variance",
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ) |>
      suppressMessages(),
    paste("The object array_manifest does not have the required columns: chr,",
          "pos, strand . Otherwise, provide a string with one of the supported",
          "human microarrays ('IlluminaHumanMethylation450k',",
          "'IlluminaHumanMethylationEPICv1', or 'IlluminaHumanMethylationEPICv2')."),
    fixed = TRUE
  )
  #### cor_threshold ####
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = "1000",
      var_method = "variance",
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ) |>
      suppressMessages(),
    "Please make sure the input cor_threshold belongs to the numeric class."
  )
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 1000,
      var_method = "variance",
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ) |>
      suppressMessages(),
    "'cor_threshold' must be a value between 0 and 1 (inclusive)",
    fixed = TRUE
  )
  #### var_method ####
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 0,
      var_method = 1,
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ) |>
      suppressMessages(),
    "Please make sure the input var_method belongs to the character class."
  )
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 0,
      var_method = c("variance", "mad"),
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ) |>
      suppressMessages(),
    "Please make sure the input var_method is a character object of length 1"
  )
  #### var_distribution ####
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 0,
      var_method = "variance",
      var_distribution = 1,
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ) |>
      suppressMessages(),
    "Please make sure the input var_distribution belongs to the character class."
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
    ) |>
      suppressMessages(),
    "Please make sure the input var_distribution is one of the following options: ultrastable, all"
  )
  #### var_threshold_percentile ####
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 0,
      var_method = "variance",
      var_distribution = "ultrastable",
      var_threshold_percentile = "0.99",
      max_distance = 1000
    ) |>
      suppressMessages(),
    "Please make sure the input var_threshold_percentile belongs to the numeric class."
  )
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 0,
      var_method = "variance",
      var_distribution = "ultrastable",
      var_threshold_percentile = 100,
      max_distance = 1000
    ) |>
      suppressMessages(),
    "'var_threshold_percentile' must be a value between 0 and 1 (inclusive)",
    fixed = TRUE
  )
  #### max_distance ####
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 0,
      var_method = "variance",
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = "1000"
    ) |>
      suppressMessages(),
    "Please make sure the input max_distance belongs to the numeric class."
  )
})

test_that("findVML works with EPICv2 probes", {
  # Set the parallel backend to use 2 workers
  doParallel::registerDoParallel(2)
  epic2_methylation_data <- RAMEN::test_methylation_data
  rownames(epic2_methylation_data) <- data.frame(IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations) |>
    dplyr::filter(chr == "chr21") |>
    dplyr::arrange(chr, pos) |>
    dplyr::slice_head(n = nrow(RAMEN::test_methylation_data)) |>
    rownames()

  VML_epic2 <- RAMEN::findVML(
    methylation_data = epic2_methylation_data,
    array_manifest = "IlluminaHumanMethylationEPICv2",
    cor_threshold = 0,
    var_method = "variance",
    var_distribution = "ultrastable",
    var_threshold_percentile = 0.99,
    max_distance = 1000
  ) |>
    suppressMessages()
  expect_true(is.list(VML_epic2))
  expect_true(is(VML_epic2$VML, "GRanges"))
  expect_equal(ncol(mcols(VML_epic2$VML)), 5)
})

test_that("findVML works with var_distribution = 'all' and mad score", {
  # Set the parallel backend to use 2 workers
  doParallel::registerDoParallel(2)
  VML_allvar <- RAMEN::findVML(
    methylation_data = RAMEN::test_methylation_data,
    array_manifest = "IlluminaHumanMethylationEPICv1",
    cor_threshold = 0,
    var_method = "mad",
    var_distribution = "all",
    var_threshold_percentile = 0.9,
    max_distance = 1000
  ) |>
    suppressMessages()
  expect_true(is.list(VML_allvar))
  expect_true(is(VML_allvar$VML, "GRanges"))
  expect_equal(ncol(mcols(VML_allvar$VML)), 5)
})

test_that("sVMPs have no correlation", {
  sVMPs <- VML_test$VML[VML_test$VML$type == "sVMP"]
  expect_true(all(is.na(sVMPs$median_correlation)))
})

test_that("correlation is computed correctly", {
  vmr_gr <- VML_test$VML[VML_test$VML$type == "VMR"]
  # Order by n_VMPs ascending, then take the last one (highest n_VMPs)
  vmr_gr <- vmr_gr[order(vmr_gr$n_VMPs)]
  VML_test_cor <- vmr_gr[length(vmr_gr)]
  # get the probes
  probes <- unlist(VML_test_cor$probes)
  methylation_subset <- RAMEN::test_methylation_data[probes, ]
  cor_matrix <- cor(t(methylation_subset))
  cor_values <- cor_matrix[lower.tri(cor_matrix)]
  median_cor <- median(cor_values)
  expect_equal(VML_test_cor$median_correlation, median_cor)
})
