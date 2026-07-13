test_that("summarizeVML output structure is correct", {
  expect_true(is.matrix(summarized_methyl_VML_test))
  expect_equal(ncol(summarized_methyl_VML_test), length(VML_test$VML))
  expect_equal(nrow(summarized_methyl_VML_test), ncol(test_methylation_data))
})

test_that("summarizeVML adds VML_index when not present", {
  # Set the parallel backend to use 2 workers
  doParallel::registerDoParallel(2)
  # Delete IDs
  VML_no_IDs <- VML_test$VML[1:5, ]
  VML_no_IDs$VML_index <- NULL
  # Run function
  summarized_no_index <- RAMEN::summarizeVML(
    VML = VML_no_IDs,
    methylation_data = RAMEN::test_methylation_data
  )
  expect_true(is.matrix(summarized_no_index))
  # All VMLs are included
  expect_true(all(colnames(summarized_no_index) %in% paste0("VML",
                                                           seq_len(length(VML_no_IDs)))
  )
  # All individuals are included
  expect_equal(nrow(summarized_no_index), ncol(RAMEN::test_methylation_data))
})

test_that("summarizeVML values are correct", {
  # Set the parallel backend to use 2 workers
  doParallel::registerDoParallel(2)
  # First for sVMPs: the summarized value should be equal to the methylation
  VML_gr <- VML_test$VML
  VML_test_summ <- VML_gr[VML_gr$type == "sVMP"][1]  # Get the first sVMP
  probe <- unlist(VML_test_summ$probes)
  expected <- RAMEN::test_methylation_data[probe, ] |>
    unlist()
  observed <- summarized_methyl_VML_test[, VML_test_summ$VML_index]
  names(observed) <- rownames(summarized_methyl_VML_test)
  expect_equal(observed, expected)
  # now for VMRs: the summarized value should be the median across probes
  VML_test_summ <- VML_gr[VML_gr$type == "VMR"][1]  # Get the first VMR
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
      VML = "a",
      methylation_data = test_methylation_data
    ),
    "Please make sure the input VML belongs to the GRanges class."
  )
  expect_error(
    RAMEN::summarizeVML(
      VML = VML_test$VML,
      methylation_data = "a"
    ),
    "Please make sure the input methylation_data belongs to the data.frame class."
  )
})
