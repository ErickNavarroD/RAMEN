#' Simulate a delta R squared null distribution of G and E effects on DNAme
#' variability
#'
#' This function simulates the delta R squared distribution under the null
#' hypothesis of G and E having no association with DNA methylation (DNAme)
#' variability through a permutation analysis. To do so, this function shuffles
#' the G and E variables in the dataset, which is followed by a the variable
#' selection and modelling steps with *selectVariables()* and *lmGE()*.These
#' steps are repeated several times as indicated in the *permutations* parameter.
#'  By using shuffled G and E data, we simulate the increase of R2 that would be
#'   observed in random data using the RAMEN methodology.
#'
#' The core pipeline from the RAMEN package identifies the best explanatory
#' model per VML. However, despite these models being winners in comparison to
#' models including any other G/E variable(s) in the dataset, some winning
#' models might perform no better than what we would expect by chance. Therefore,
#' the goal of this function is to create a distribution of increase in R2
#' under the null hypothesis of G and E having no associations with DNAme.
#' The null distribution is obtained through shuffling the G and E variables
#' in a given dataset and conducting the variable selection and G/E model
#' selection. That way, we can simulate how much additional variance would be
#' explained by the models defined as winners by the RAMEN methodology in a
#' scenario where the G and E associations with DNAme are randomized. This
#' distribution can be then used to filter out winning models in the
#' non-shuffled dataset that do not add more to the explained variance of the
#' basal model than what randomized data do.
#'
#' Under the assumption that after adjusting for the concomitant variables all
#' VML across the genome follow the same behavior regarding an increment of
#' explained variance with randomized G and E data, we can pool the delta R
#' squared values from all VML to create a null distribution taking advantage
#' of the high number of VML in the dataset. This assumption decreases
#' significantly the number of permutations required to create a null
#' distribution and reduces the computational time. For further information
#' please read the RAMEN paper (https://doi.org/10.1186/s13059-025-03864-4).
#'
#' @inheritParams selectVariables
#' @inheritParams lmGE
#' @param permutations Numer of permutation analyses to run.
#'
#' @return A data frame with the following columns:
#'  - VML_index: The unique ID of the VML.
#'  - model_group: The group to which the winning model belongs to (i.e., G, E,
#'  G+E or GxE)
#'  - tot_r_squared: R squared of the winning model
#'  - R2_difference: the increase in R squared obtained by including the G/E
#'  variable(s) from the winning model (i.e., the R squared difference between
#'  the winning model and the model only with the concomitant variables
#'  specified in *covariates*; tot_r_squared - basal_rsquared in the lmGE output)
#'  - AIC_difference/BIC_difference: the AIC/BIC difference between the winning
#'  model and the model only with the concomitant variables specified in
#'  *covariates*; BIC/AIC - basal_BIC/basal_BIC in the lmGE output)
#'
#' @importFrom foreach %do%
#' @export
#' @examples
#' # Evaluate sequentially
#' foreach::registerDoSEQ()
#' ## Find VML in test data
#' VML <- findVML(
#'   methylation_data = test_methylation_data,
#'   array_manifest = "IlluminaHumanMethylationEPICv1",
#'   cor_threshold = 0,
#'   var_method = "variance",
#'   var_distribution = "ultrastable",
#'   var_threshold_percentile = 0.99,
#'   max_distance = 1000
#' )
#' ## Find cis SNPs around VML
#' VML_with_cis_snps <- findCisSNPs(
#'   # Use only 5 for demonstration purposes
#'   VML = VML$VML[1:5, ],
#'   genotype_information = test_genotype_information,
#'   distance = 1e6
#' )
#'
#' ## Summarize methylation levels in VML
#' summarized_methyl_VML <- summarizeVML(
#'   methylation_data = test_methylation_data,
#'   VML = VML_with_cis_snps
#' )
#'
#' ## Simulate null distribution of G and E contributions on DNAme variability
#' ## We will only run one permutation for demonstration purposes
#' null_dist <- nullDistGE(
#'   VML_wSNPs = VML_with_cis_snps,
#'   genotype_matrix = test_genotype_matrix,
#'   environmental_matrix = test_environmental_matrix,
#'   summarized_methyl_VML = summarized_methyl_VML,
#'   # Use one permutation for demonstration purposes
#'   permutations = 1,
#'   covariates = test_covariates,
#'   seed = 1,
#'   model_selection = "AIC"
#' )
nullDistGE <- function(VML_wSNPs,
                       genotype_matrix,
                       environmental_matrix,
                       summarized_methyl_VML,
                       permutations = 5,
                       covariates = NULL,
                       seed = NULL,
                       model_selection = "AIC") {
  #### Argument check ####
  ## Check that genotype_matrix, environmental_matrix, and covariates (in case
  ## it is provided) have only numeric values and no NA, NaN, Inf values
  argument_check(VML_wSNPs, "GRanges")
  columns_exist(S4Vectors::mcols(VML_wSNPs), c("VML_index", "SNP"))
  argument_check(VML_wSNPs$SNP, "list")
  argument_check(VML_wSNPs$VML_index, "character")
  argument_check(genotype_matrix, "matrix")
  argument_check(environmental_matrix, "matrix")
  argument_check(summarized_methyl_VML, "matrix")
  argument_check(permutations, "numeric")
  if (permutations < 1) stop("Please provide a permutation number >= 1")
  if (!is.null(covariates)) argument_check(covariates, "matrix")
  if (!is.null(seed)) argument_check(seed, "numeric")
  argument_char_options(object = model_selection, options = c("AIC", "BIC"))
  # Matrices have only finite numeric values
  finite_numeric_check(genotype_matrix)
  finite_numeric_check(environmental_matrix)
  finite_numeric_check(summarized_methyl_VML)
  if (!is.null(covariates)) finite_numeric_check(covariates)

  #### Shuffle data ####
  # Set shuffle order
  if (!is.null(seed)) set.seed(seed)
  # Initialize permutation object
  permutation_order <- data.frame(sample(rownames(summarized_methyl_VML),
                                         size = length(rownames(summarized_methyl_VML))
  ))
  if (permutations > 1) {
    # Append order of other permutations
    for (i in 1:(permutations - 1)) {
      permutation_order <- cbind(
        permutation_order,
        data.frame(sample(rownames(summarized_methyl_VML),
                          size = length(rownames(summarized_methyl_VML))
        )))
    }
  }
  colnames(permutation_order) <- 1:permutations

  # Put the environmental and genotype matrix in the same order to the
  # summarized VML object
  genotype_matrix <- genotype_matrix[, rownames(summarized_methyl_VML)]
  environmental_matrix <- environmental_matrix[rownames(summarized_methyl_VML), ]

  #### Run permutation ####
  null_dist <- foreach::foreach(i = 1:permutations, .combine = rbind) %do% {
    message("Starting permutation ", i, " of ", permutations)
    # Change order of samples
    permutated_genotype <- genotype_matrix[, permutation_order[, i],
                                           drop = FALSE]
    permutated_environment <- environmental_matrix[permutation_order[, i], ,
                                                   drop = FALSE]
    # Assign previous row and colnames - break relationships, i.e., shuffle
    colnames(permutated_genotype) <- colnames(genotype_matrix)
    rownames(permutated_environment) <- rownames(environmental_matrix)

    # Run RAMEN
    message("Starting variable selection of permutation ", i, " of ", permutations)
    selected_variables <- RAMEN::selectVariables(
      VML_wSNPs = VML_wSNPs,
      genotype_matrix = permutated_genotype,
      environmental_matrix = permutated_environment,
      covariates = covariates,
      summarized_methyl_VML = summarized_methyl_VML,
      seed = 1
    )

    message("Starting lmGE in permutation ", i, " of ", permutations)
    lmGE_res <- RAMEN::lmGE(
      selected_variables = selected_variables,
      summarized_methyl_VML = summarized_methyl_VML,
      genotype_matrix = permutated_genotype,
      environmental_matrix = permutated_environment,
      covariates = covariates,
      model_selection = model_selection
    )
    if (model_selection == "AIC") {
      results_perm <- data.frame(
        VML_index = lmGE_res$VML_index,
        model_group = lmGE_res$model_group,
        tot_r_squared = lmGE_res$tot_r_squared,
        R2_difference = lmGE_res$tot_r_squared - lmGE_res$basal_rsquared,
        AIC_difference = lmGE_res$AIC - lmGE_res$basal_AIC
      )
    } else if (model_selection == "BIC") {
      results_perm <- data.frame(
        VML_index = lmGE_res$VML_index,
        model_group = lmGE_res$model_group,
        tot_r_squared = lmGE_res$tot_r_squared,
        R2_difference = lmGE_res$tot_r_squared - lmGE_res$basal_rsquared,
        BIC_difference = lmGE_res$BIC - lmGE_res$basal_BIC
      )
    }
    message("Wrapping up permutation ", i, " of ", permutations)
    results_perm$permutation <- i # add the number of permutation
    results_perm
  }
  #### Return results ####
  return(null_dist)
}
