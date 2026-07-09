test_that("findCisSNPs output structure is correct", {
  expect_true(is(VML_cis_snps_test, "GRanges"))
  expect_equal(ncol(S4Vectors::mcols(VML_cis_snps_test)),
               ncol(S4Vectors::mcols(VML_test$VML)) + 2)
  expect_equal(length(VML_cis_snps_test), length(VML_test$VML))
  expect_true(all(
    c(colnames(S4Vectors::mcols(VML_test$VML)),
      "surrounding_SNPs",
      "SNP") %in% colnames(S4Vectors::mcols(VML_cis_snps_test))
  ))
})

test_that("findCisSNPs adds a VML index when it is not present", {
  # Set the parallel backend to use 2 workers
  doParallel::registerDoParallel(2)
  # Delete IDs
  VML_no_IDs <- VML_test$VML[1:5, ]
  VML_no_IDs$VML_index <- NULL
  # Run function
  VML_cis_snps_noID <- RAMEN::findCisSNPs(
    VML = VML_no_IDs,
    genotype_information = RAMEN::test_genotype_information,
    distance = 1e+06
  ) |>
    suppressMessages()
  expect_true("VML_index" %in% colnames(S4Vectors::mcols(VML_cis_snps_noID)))
})

test_that("findCisSNPs throws errors when expected", {
  expect_error(
    RAMEN::findCisSNPs(
      VML = data.frame(VML_test$VML),
      genotype_information = RAMEN::test_genotype_information,
      distance = 1e+06
    ) |>
      suppressMessages(),
    "Please make sure the input VML belongs to the GRanges class",
    fixed = TRUE
  )
  expect_error(
    RAMEN::findCisSNPs(
      VML = VML_test$VML,
      genotype_information = RAMEN::test_genotype_information |>
        dplyr::select(-CHROM),
      distance = 1e+06
    ) |>
      suppressMessages(),
    "The object genotype_information does not have the required columns: CHROM, POS, ID .",
    fixed = TRUE
  )
  expect_error(
    RAMEN::findCisSNPs(
      VML = VML_test$VML,
      genotype_information = "a",
      distance = 1e+06
    ) |>
      suppressMessages(),
    "Please make sure the input genotype_information belongs to the data.frame class."
  )
  expect_error(
    RAMEN::findCisSNPs(
      VML = VML_test$VML,
      genotype_information = RAMEN::test_genotype_information,
      distance = "a"
    ) |>
      suppressMessages(),
    "Please make sure the input distance belongs to the numeric class."
  )
})

test_that("findCisSNPs returns the right number of cis SNPs", {
  # Set the parallel backend to use 2 workers
  doParallel::registerDoParallel(2)
  VML_vanilla <- GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1000, end = 2000),
    VML_index = "1",
    type = "VMR"
  )
  genot_info_test <- data.frame(
    CHROM = c("chr1", "chr1", "chr1", "chr1"),
    POS = c(1, 500, 2500, 4000),
    ID = c("rs1", "rs2", "rs3", "rs4")
  )
  test_1 <- RAMEN::findCisSNPs(
    VML = VML_vanilla,
    genotype_information = genot_info_test,
    distance = 1
  ) |>
    suppressMessages()
  test_500 <- RAMEN::findCisSNPs(
    VML = VML_vanilla,
    genotype_information = genot_info_test,
    distance = 500
  ) |>
    suppressMessages()
  test_1000 <- RAMEN::findCisSNPs(
    VML = VML_vanilla,
    genotype_information = genot_info_test,
    distance = 1000
  ) |>
    suppressMessages()
  test_2000 <- RAMEN::findCisSNPs(
    VML = VML_vanilla,
    genotype_information = genot_info_test,
    distance = 2000
  ) |>
    suppressMessages()
  expect_equal(test_1$surrounding_SNPs, 0) # no SNPs are within 1 bp
  expect_equal(test_500$surrounding_SNPs, 2) # only rs2 and rs3 are within 500bp
  expect_equal(test_1000$surrounding_SNPs, 3) # rs1, rs2 and rs3 are within 1000bp
  expect_equal(test_2000$surrounding_SNPs, 4) # all 4 snps are within 2000bp
})
