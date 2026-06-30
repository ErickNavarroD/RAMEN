test_that("findCisSNPs output structure is correct", {
  expect_true(is.data.frame(VML_cis_snps_test))
  expect_equal(ncol(VML_cis_snps_test), ncol(VML_test$VML) + 2)
  expect_equal(nrow(VML_cis_snps_test), nrow(VML_test$VML))
  expect_true(all(
    c(colnames(VML_test$VML), "surrounding_SNPs", "SNP") %in%
      colnames(VML_cis_snps_test)
  ))
})

test_that("findCisSNPs adds a VML index when it is not present", {
  VML_cis_snps_noID <- RAMEN::findCisSNPs(
    VML_df = VML_test$VML |>
      dplyr::select(-VML_index),
    genotype_information = RAMEN::test_genotype_information,
    distance = 1e+06
  ) |>
    suppressMessages()
  expect_true("VML_index" %in% colnames(VML_cis_snps_noID))
})

test_that("findCisSNPs throws errors when expected", {
  expect_error(
    RAMEN::findCisSNPs(
      VML_df = VML_test$VML |>
        dplyr::select(-seqnames),
      genotype_information = RAMEN::test_genotype_information,
      distance = 1e+06
    ) |>
      suppressMessages(),
    "Please make sure the VML_df object has the required columns with the appropiate names (check documentation for further information)",
    fixed = TRUE
  )
  expect_error(
    RAMEN::findCisSNPs(
      VML_df = VML_test$VML,
      genotype_information = RAMEN::test_genotype_information |>
        dplyr::select(-CHROM),
      distance = 1e+06
    ) |>
      suppressMessages(),
    "Please make sure the genotype_information object has the required columns with the appropiate names (check documentation for further information)",
    fixed = TRUE
  )
  expect_error(
    RAMEN::findCisSNPs(
      VML_df = "a",
      genotype_information = RAMEN::test_genotype_information,
      distance = 1e+06
    ) |>
      suppressMessages(),
    "Please make sure the VML_df object is a data frame.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::findCisSNPs(
      VML_df = VML_test$VML,
      genotype_information = "a",
      distance = 1e+06
    ) |>
      suppressMessages(),
    "Please make sure the genotype_information object is a data frame."
  )
})

test_that("findCisSNPs returns the right number of cis SNPs", {
  VML_vanilla <- data.frame(
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
    VML_df = VML_vanilla,
    genotype_information = genot_info_test,
    distance = 1
  ) |>
    suppressMessages()
  test_500 <- RAMEN::findCisSNPs(
    VML_df = VML_vanilla,
    genotype_information = genot_info_test,
    distance = 500
  ) |>
    suppressMessages()
  test_1000 <- RAMEN::findCisSNPs(
    VML_df = VML_vanilla,
    genotype_information = genot_info_test,
    distance = 1000
  ) |>
    suppressMessages()
  test_2000 <- RAMEN::findCisSNPs(
    VML_df = VML_vanilla,
    genotype_information = genot_info_test,
    distance = 2000
  ) |>
    suppressMessages()
  expect_equal(test_1$surrounding_SNPs, 0) # no SNPs are within 1 bp
  expect_equal(test_500$surrounding_SNPs, 2) # only rs2 and rs3 are within 500bp
  expect_equal(test_1000$surrounding_SNPs, 3) # rs1, rs2 and rs3 are within 1000bp
  expect_equal(test_2000$surrounding_SNPs, 4) # all 4 snps are within 2000bp
})
