#' Selection of relevant environment and genotype variables associated with Variably Methylated Loci (VML)
#'
#' For each VML, this function selects potentially relevant genotype and environmental variables associated with DNA methylation levels of said VML using LASSO. See details below for more information.
#'
#' selectVariables() uses LASSO, which is an embedded variable selection method that penalizes models that are more complex (i.e., that contain more variables) in favor of simpler models (i.e. that contain less variables), but not at the expense of reducing predictive power. Using LASSO's variable screening property (with high probability, the LASSO estimated model includes the substantial covariates and drops the redundant ones) this function selects genotype and environment variables with potential relevance in the Variable Methylated Loci (VML) dataset (see also Bühlmann and van de Geer, 2011). For each VML, LASSO is run three times: 1) including only the genotype variables for the selection step, 2) including only the environmental variables for the selection step, and 3) Including both the genotype and environmental variables in the selection step. This is done to ensure that the function captures the variables that are relevant within their own category (e.g., SNPs that are strongly associated with the DNAme levels of a VML in the presence of the rest of the SNPs) or in the presence of the variables of the other category (e.g. SNPs that are strongly associated with the DNAme levels of a VML in the presence of the rest of BOTH the SNPs AND environmental variables). Every time LASSO is run, the basal covariates (i.e., concomitant variables )indicated in the argument *covariates* are not penalized (i.e., those variables are always included in the models and their coefficients are not subjected to shrinkage). That way, only the most promising E and G variables in the presence of the concomitant variables will be selected.
#'
#' Each LASSO model uses a tuned lambda that minimizes the 5-fold cross-validation error within its corresponding data. This function uses the lambda.min value in contrast to lambda.1se because its goal within the RAMEN package is to use LASSO to reduce the number of variables that are going to be used next for fitting pairwise interaction models in *lmGE()*. Since at this step variables are being selected based only on main effects, it is preferable to cast a "wider net" and select a slightly higher number of variables that could potentially have a strong interaction effect when paired with another variable. Furthermore, since in this case LASSO is being used as a screening procedure to select variables that will be fit separately in independent models and compared, the overfitting issue of using lambda.min does not impose a big concern. After finding the best lambda value, the sequence of models is fit by coordinate descent using *glmnet()*. Random numbers in this function are created during the lambda cross validation and the LASSO stages. Setting a seed is highly encouraged for result reproducibility using the *seed* argument. Please note that setting a seed inside of this function modifies the seed globally (which is R's default behavior).
#'
#' #' This function supports parallel computing for increased speed. To do so, you have to set the parallel back-end
#' in your R session before running the function  (e.g., *doParallel::registerDoParallel(4)*). After that, the function can be run as usual. It is recommended to also set options(future.globals.maxSize= +Inf). Please make sure that your data has no NAs and it's all numerical, since the LASSO implementation we use does not support missing or non-numerical values.
#'
#' Note: If you want to conduct the variable selection step only in one data set (i.e., only in the genotype), you can set the argument *environmental_matrix = NULL*.
#'
#'
#' @param VML_df A data frame converted from a GRanges object. Recommended to use the output of *RAMEN::findCisSNPs()*. Must have one VML per row, and contain the following columns: "VML_index" (a unique ID for each VML in VML_df AS CHARACTERS) and "SNP" (a column with a list as observation, containing the name of the SNPs surrounding the corresponding VML).  The SNPs contained in the "SNP" column must be present in the object that is indicated in the genotype_matrix argument, and it must contain all the VML contained in summarized_methyl_VML. VML with no surrounding SNPs must have an empty list in the SNP column (either list(NULL), list(NA), list("") or list(character(0)) ).
#' @param environmental_matrix A matrix of environmental variables. Only numeric values are supported. In case of factor variables, it is recommended to encode them as numbers or re-code them into dummy variables if there are more than two levels. Columns must correspond to environmental variables and rows to individuals. Row names must be the individual IDs.
#' @param genotype_matrix A matrix of number-encoded genotypes. Columns must correspond to samples, and rows to SNPs. We suggest using a gene-dosage model, which would encode the SNPs ordinally depending on the genotype allele charge, such as 2 (AA), 1 (AB) and 0 (BB). The column names must correspond with individual IDs.
#' @param summarized_methyl_VML A data frame containing each individual's VML summarized methylation. It is suggested to use the output of RAMEN::summarizeVML().Rows must reflects individuals, and columns VML The names of the columns must correspond to the index of said VML, and it must match the index of VML_df$VML_index. The names of the rows must correspond to the sample IDs, and must match with the IDs of the other matrices.
#' @param covariates A matrix containing the covariates (i.e., concomitant variables / variables that are not the ones you are interested in) that will be adjusted for in the final GxE models (e.g., cell type proportions, age, etc.). Each column should correspond to a covariate and each row to an individual. Row names must correspond to the individual IDs.
#' @param seed An integer number that initializes a pseudo-random number generator. Random numbers in this function are created during the lambda cross validation and the LASSO stages. Setting a seed is highly encouraged for result reproducibility. **Please note that setting a seed in this function modifies the seed globally**.
#'
#' @return A data frame with three columns:
#'  - VML_index: Unique VML ID.
#'  - selected_genot: Column containing lists as values with the selected SNPs.
#'  - selected_env: Column containing lists as values with the selected environmental variables.
#'
#' @importFrom doRNG %dorng%
#' @export
#'
#' @examples
#' ## Find VML in test data
#' VML <- RAMEN::findVML(
#'    methylation_data = RAMEN::test_methylation_data,
#'    array_manifest = "IlluminaHumanMethylationEPICv1",
#'    cor_threshold = 0,
#'    var_method = "variance",
#'    var_distribution = "ultrastable",
#'    var_threshold_percentile = 0.99,
#'    max_distance = 1000
#'    )
#' ## Find cis SNPs around VML
#' VML_with_cis_snps <- RAMEN::findCisSNPs(
#'   VML_df = VML$VML,
#'   genotype_information = RAMEN::test_genotype_information,
#'   distance = 1e6
#'   )
#'
#' ## Summarize methylation levels in VML
#' summarized_methyl_VML <- RAMEN::summarizeVML(
#'  methylation_data = RAMEN::test_methylation_data,
#'  VML_df = VML_with_cis_snps
#'  )
#'
#'  ## Select relevant genotype and environmental variables
#'  selected_vars <- RAMEN::selectVariables(
#'    VML_df = VML_with_cis_snps,
#'    genotype_matrix = RAMEN::test_genotype_matrix,
#'    environmental_matrix = RAMEN::test_environmental_matrix,
#'    covariates = RAMEN::test_covariates,
#'    summarized_methyl_VML = summarized_methyl_VML,
#'    seed = 1
#'  )
#'
selectVariables <- function(VML_df,
                            genotype_matrix,
                            environmental_matrix,
                            covariates = NULL,
                            summarized_methyl_VML,
                            seed = NULL) {
  # Arguments check
  ## Check that genotype_matrix, environmental_matrix, covariate matrix (in case it is provided) and summarized_methyl_VML have the same samples
  if (!all(rownames(summarized_methyl_VML) %in% colnames(genotype_matrix))) stop("Individual IDs in summarized_methyl_VML do not match individual IDs in genotype_matrix")
  if (!all(rownames(summarized_methyl_VML) %in% rownames(environmental_matrix))) stop("Individual IDs in summarized_methyl_VML do not match individual IDs in environmental_matrix")
  if (!is.null(covariates)) {
    if (!all(rownames(summarized_methyl_VML) %in% rownames(covariates))) stop("Individual IDs in summarized_methyl_VML do not match individual IDs in the covariates matrix")
  }
  ## Check that VML_df has index and SNP column
  if (!all(c("VML_index", "SNP") %in% colnames(VML_df))) stop("Please make sure the VML data frame (VML_df) contains the columns 'SNP' and 'VML_index'.")
  ## Check that the SNP column on VML_df is a list
  if (!is.list(VML_df$SNP)) stop("Please make sure the 'SNP' column in VML_df is a column containing lists as values")
  if (!is.character(VML_df$VML_index)) stop("Please make sure the 'VML_index' column in VML_df is a column of characters")
  ## Check that genotype, environment and covariates are matrices
  if (!is.matrix(genotype_matrix)) stop("Please make sure the genotype data is provided as a matrix.")
  if (!is.null(environmental_matrix)) {
    if (!is.matrix(environmental_matrix)) stop("Please make sure the environmental data is provided as a matrix.")
  }
  if (!is.null(covariates)) {
    if (!is.matrix(covariates)) stop("Please make sure the covariates data is provided as a matrix.")
  }
  ## Check that genotype_matrix, environmental_matrix, and covariates (in case
  ## it is provided) have only numeric values and no NA, NaN, Inf values
  if (
    sum(sapply(genotype_matrix, is.na)) > 0 ||
    sum(sapply(genotype_matrix, is.nan)) > 0 ||
    sum(!sapply(genotype_matrix, is.numeric)) > 0 ||
    sum(sapply(genotype_matrix, is.infinite)) > 0
    ) stop (
      "Please make sure the genotype matrix contains only finite numeric values."
      )
  if (
    sum(sapply(environmental_matrix, is.na)) > 0 ||
    sum(sapply(environmental_matrix, is.nan)) > 0 ||
    sum(!sapply(environmental_matrix, is.numeric)) > 0 ||
    sum(sapply(environmental_matrix, is.infinite)) > 0
  ) stop (
    "Please make sure the environmental matrix contains only finite numeric values."
  )
  if (!is.null(covariates)) {
    if (
      sum(sapply(covariates, is.na)) > 0 ||
      sum(sapply(covariates, is.nan)) > 0 ||
      sum(!sapply(covariates, is.numeric)) > 0 ||
      sum(sapply(covariates, is.infinite)) > 0
    ) stop (
      "Please make sure the covariates matrix contains only finite numeric values."
    )
  }

  if (
    sum(sapply(summarized_methyl_VML, is.na)) > 0 ||
    sum(sapply(summarized_methyl_VML, is.nan)) > 0 ||
    sum(!sapply(summarized_methyl_VML, is.numeric)) > 0 ||
    sum(sapply(summarized_methyl_VML, is.infinite)) > 0
  ) stop (
    "Please make sure the summarized_methyl_VML matrix or data frame contains only finite numeric values."
  )

  ## Set the seed
  if (!is.null(seed)) set.seed(seed)

  lasso_results <- foreach::foreach(VML_i = iterators::iter(VML_df, by = "row"), .combine = "rbind") %dorng% {
    # Select summarized VML information
    summVMLi <- summarized_methyl_VML %>%
      dplyr::select(VML_i$VML_index)
    ## Prepare data
    # subset the genotyping data and match genotype, environment and DNAme IDs
    if (VML_i$SNP %in% list(NULL) | # Catch VML with no surrounding SNPs
        VML_i$SNP %in% list("") |
        VML_i$SNP %in% list(NA) |
        VML_i$SNP %in% list(character(0))) {
      genot_VMLi <- c()
      any_snp <- FALSE
    } else if (length(VML_i$SNP[[1]]) == 1) { # Special case of sub-setting if SNP is only one because the result is a vector and not a matrix
      genot_VMLi <- genotype_matrix[unlist(VML_i$SNP), rownames(summVMLi)] %>%
        as.matrix()
      colnames(genot_VMLi) <- VML_i$SNP[[1]]
      any_snp <- TRUE
    } else {
      genot_VMLi <- genotype_matrix[unlist(VML_i$SNP), rownames(summVMLi)] %>%
        t()
      any_snp <- TRUE
    }
    if (ncol(environmental_matrix) == 1) {
      environ_VMLi <- environmental_matrix[rownames(summVMLi), ] %>%
        as.matrix()
      colnames(environ_VMLi) <- colnames(environmental_matrix)
    } else {
      environ_VMLi <- environmental_matrix[rownames(summVMLi), ]
    }
    environ_genot_VMLi <- cbind(genot_VMLi, environ_VMLi)
    # Bind covariates data
    if (!is.null(covariates)) {
      if (ncol(covariates) == 1) {
        covariates_VMLi <- covariates[rownames(summVMLi), ] %>% # Match the covariates dataset with the VML information
          as.matrix()
        colnames(covariates_VMLi) <- colnames(covariates)
      } else {
        covariates_VMLi <- covariates[rownames(summVMLi), ]
      }
      genot_VMLi <- cbind(genot_VMLi, covariates_VMLi)
      environ_VMLi <- cbind(environ_VMLi, covariates_VMLi)
      environ_genot_VMLi <- cbind(environ_genot_VMLi, covariates_VMLi)
      ncol_covariates <- ncol(covariates_VMLi)
    } else {
      ncol_covariates <- 0
    }

    ### Run LASSOs
    ## Genotype only
    # Get coefficients with the optimal lambda found by k-fold cross-validation
    if (any_snp) { # Catch cases when VML dont have surrounding genotyped SNPs
      coef_genot <- stats::coef(
        glmnet::cv.glmnet(
          x = genot_VMLi, # Variables
          y = summVMLi[, VML_i$VML_index], # Response
          alpha = 1,
          nfolds = 5,
          penalty.factor = c(
            rep(1, ncol(genot_VMLi) - ncol_covariates),
            rep(0, ncol_covariates)
          )
        ), # Unpenalize the variables in covariates (i.e., force LASSO to keep them in all the situations)
        s = "lambda.min"
      )
      # Select the variables with a coefficient > 0
      coef_genot <- coef_genot[abs(coef_genot[, 1]) > 0, ]
      selected_vars_genot <- names(coef_genot)[-1]
      selected_vars_genot <- selected_vars_genot[!selected_vars_genot %in% colnames(covariates)] # Remove covariates from selected variables
    } else {
      selected_vars_genot <- character(0)
    }

    # Environment only
    # Get coefficients with the optimal lambda found by k-fold cross-validation
    if (!is.null(environ_VMLi)) { # catch scenario where users would not add environmental variables
      coef_env <- stats::coef(
        glmnet::cv.glmnet(
          x = environ_VMLi, # Variables
          y = summVMLi[, VML_i$VML_index], # Response
          alpha = 1,
          nfolds = 5,
          penalty.factor = c(
            rep(1, ncol(environ_VMLi) - ncol_covariates), # Unpenalize the variables in covariates (i.e., force LASSO to keep them in all the situations)
            rep(0, ncol_covariates)
          )
        ),
        s = "lambda.min"
      )
      # Select the variables with a coefficient > 0
      coef_env <- coef_env[abs(coef_env[, 1]) > 0, ]
      selected_vars_env <- names(coef_env)[-1] # Remove the intercept from the variables
      selected_vars_env <- selected_vars_env[!selected_vars_env %in% colnames(covariates)] # Remove covariates from selected variables

      if (any_snp) {
        # Joint (environment + genotype) only when we have Genotype and Environmental variables.
        # Get coefficients with the optimal lambda found by k-fold cross-validation
        coef_joint <- stats::coef(
          glmnet::cv.glmnet(
            x = environ_genot_VMLi, # Variables
            y = summVMLi[, VML_i$VML_index], # Response
            alpha = 1,
            nfolds = 5,
            penalty.factor = c(
              rep(1, ncol(environ_genot_VMLi) - ncol_covariates),
              rep(0, ncol_covariates)
            )
          ), # Unpenalize the variables in covariates (i.e., force LASSO to keep them in all the situations)
          s = "lambda.min"
        )
        # Select the variables with an abs(coefficient) > 0
        coef_joint <- coef_joint[abs(coef_joint[, 1]) > 0, ]
        selected_vars_joint <- names(coef_joint)[-1] # Remove the intercept from the variables
        selected_vars_joint <- selected_vars_joint[!selected_vars_joint %in% colnames(covariates)] # Remove covariates from selected variables
      } else {
        selected_vars_joint <- character(0)
      }
    } else {
      selected_vars_env <- character(0)
      selected_vars_joint <- character(0)
    }

    # Merge results
    selected_union_genot <- c(selected_vars_genot, selected_vars_joint) %>%
      unique() %>%
      dplyr::setdiff(colnames(environ_VMLi)) # Remove environmental variables and covariates from the joint selection
    selected_union_env <- c(selected_vars_env, selected_vars_joint) %>%
      unique() %>%
      dplyr::setdiff(colnames(genot_VMLi)) # Remove genotype variables and covariates from the joint selection

    ### Create final data frame
    selected_variables_final <- data.frame(
      VML_index = VML_i$VML_index
    )
    selected_variables_final$selected_genot <- list(selected_union_genot)
    selected_variables_final$selected_env <- list(selected_union_env)
    selected_variables_final
  }

  return(data.frame(lasso_results))
}
