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
  expect_true(is.data.frame(VML_test$VML))
  expect_true(is.data.frame(VML_test$highly_variable_probes))
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
  expect_error(
    RAMEN::findVML(
      methylation_data = RAMEN::test_methylation_data,
      array_manifest = "IlluminaHumanMethylationEPICv1",
      cor_threshold = 1000,
      var_method = "ultrastable",
      var_distribution = "ultrastable",
      var_threshold_percentile = 0.99,
      max_distance = 1000
    ) |>
      suppressMessages(),
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
    ) |>
      suppressMessages(),
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
    ) |>
      suppressMessages(),
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
    ) |>
      suppressMessages()
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
    ) |>
      suppressMessages(),
    "The method must be either 'mad' or 'variance'. Please select one of those options"
  )
})

test_that("findVML works with EPICv2 probes", {
  epic2_methylation_data <- RAMEN::test_methylation_data
  rownames(epic2_methylation_data) <- data.frame(IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations) |>
    dplyr::filter(chr == "chr21") |>
    dplyr::arrange(chr, pos) |> # Make sure to extract neighbouring probes to have VML
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
  expect_true(is.data.frame(VML_epic2$VML))
  expect_equal(ncol(VML_epic2$VML), 10)
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
  ) |>
    suppressMessages()
  expect_true(is.list(VML_allvar))
  expect_true(is.data.frame(VML_allvar$VML))
  expect_equal(ncol(VML_allvar$VML), 10)
})

test_that("sVMPs have no correlation", {
  sVMPs <- VML_test$VML |>
    dplyr::filter(type == "sVMP") |>
    dplyr::pull(median_correlation)
  expect_true(all(is.na(sVMPs)))
})

test_that("correlation is computed correctly", {
  VML_test_cor <- VML_test$VML |>
    dplyr::filter(type == "VMR") |>
    dplyr::arrange(n_VMPs) |>
    dplyr::slice_tail(n = 1) # Get the VMR with highest number of HVPs
  probes <- unlist(VML_test_cor$probes)
  methylation_subset <- RAMEN::test_methylation_data[probes, ]
  cor_matrix <- cor(t(methylation_subset))
  cor_values <- cor_matrix[lower.tri(cor_matrix)]
  median_cor <- median(cor_values)
  expect_equal(VML_test_cor$median_correlation, median_cor)
})
