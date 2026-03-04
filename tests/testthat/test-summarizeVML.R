test_that("summarizeVML output structure is correct", {
  expect_true(is.data.frame(summarized_methyl_VML_test))
  expect_equal(ncol(summarized_methyl_VML_test), nrow(VML_test$VML))
  expect_equal(nrow(summarized_methyl_VML_test), ncol(test_methylation_data))
})

test_that("summarizeVML adds VML_index when not present", {
  VML_no_index <- VML_test$VML |>
    dplyr::select(-VML_index)
  summarized_no_index <- RAMEN::summarizeVML(
    VML_df = VML_no_index,
    methylation_data = RAMEN::test_methylation_data
  )
  expect_true(is.data.frame(summarized_no_index))
  expect_true(all(
    colnames(summarized_no_index) %in% paste0("VML", seq_len(nrow(VML_no_index)))
  ))
  expect_equal(nrow(summarized_no_index), ncol(RAMEN::test_methylation_data))
})

test_that("summarizeVML values are correct", {
  # First for sVMPs: the summarized value should be equal to the methylation
  VML_test_summ <- VML_test$VML |>
    dplyr::filter(type == "sVMP") |>
    dplyr::slice_head(n = 1) # Get the first sVMP
  probe <- unlist(VML_test_summ$probes)
  expected <- RAMEN::test_methylation_data[probe, ] |> unlist()
  observed <- summarized_methyl_VML_test[, VML_test_summ$VML_index]
  names(observed) <- rownames(summarized_methyl_VML_test)
  expect_equal(observed, expected)
  # now for VMRs: the summarized value should be the median across probes
  VML_test_summ <- VML_test$VML |>
    dplyr::filter(type == "VMR") |>
    dplyr::slice_head(n = 1) # Get the first VMR
  probes <- unlist(VML_test_summ$probes)
  expected <- apply(
    RAMEN::test_methylation_data[probes, ],
    2,
    median
  )
  observed <- summarized_methyl_VML_test[, VML_test_summ$VML_index]
  names(observed) <- rownames(summarized_methyl_VML_test)
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
      VML_df = VML_test$VML,
      methylation_data = "a"
    ),
    "Please make sure the methylation data is a data frame or matrix with samples as columns and probes as rows."
  )
})

test_that("summarizeVML works when methylation_data is a matrix", {
  summarized_methyl_VML_matrix <- RAMEN::summarizeVML(
    VML_df = VML_test$VML,
    methylation_data = as.matrix(test_methylation_data)
  )
  expect_true(is.data.frame(summarized_methyl_VML_matrix))
  expect_equal(ncol(summarized_methyl_VML_matrix), nrow(VML_test$VML))
  expect_equal(nrow(summarized_methyl_VML_matrix), ncol(test_methylation_data))
  svmp <- VML_test$VML |>
    dplyr::filter(type == "sVMP") |>
    dplyr::slice_head(n = 1) # Get the first sVMP
  expect_equal(
    summarized_methyl_VML_matrix[, svmp$VML_index],
    test_methylation_data[unlist(svmp$probes), ] |> unlist() |> unname()
  )
})
