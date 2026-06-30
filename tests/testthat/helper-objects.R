suppressMessages(VML_test <- RAMEN::findVML(
  methylation_data = RAMEN::test_methylation_data,
  array_manifest = "IlluminaHumanMethylationEPICv1",
  cor_threshold = 0,
  var_method = "variance",
  var_distribution = "ultrastable",
  var_threshold_percentile = 0.99,
  max_distance = 1000
))

summarized_methyl_VML_test <- RAMEN::summarizeVML(
  VML_df = VML_test$VML,
  methylation_data = test_methylation_data
)

suppressMessages(VML_cis_snps_test <- RAMEN::findCisSNPs(
  VML_df = VML_test$VML,
  genotype_information = RAMEN::test_genotype_information,
  distance = 1e+06
))

suppressMessages(selected_variables_test <- RAMEN::selectVariables(
  VML_df = VML_cis_snps_test,
  genotype_matrix = RAMEN::test_genotype_matrix,
  environmental_matrix = RAMEN::test_environmental_matrix,
  covariates = RAMEN::test_covariates,
  summarized_methyl_VML = summarized_methyl_VML_test,
  seed = 1
))
