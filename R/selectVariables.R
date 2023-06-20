#' Selection of substantial environment and genotype variables for Variable Methylated Regions
#'
#' For each Variable Methylated Region, this functoon selects the relevant genotype and environment variables using LASSO. OPTIONAL: This function supports parallel computing for increased speed. To do so, you have to set the parallel backend in your R session BEFORE running the function (e.g., doFuture::registerDoFuture()) and then the evaluation strategy (e.g., future::plan(multisession)). After that, the function can be run as usual.
#'
#' This function uses LASSO, which is an embedded variable selection method that penalizes models that are more complex (i.e., that contain more variables) in favor of simpler models (i.e. that contain less variables), but not at the expense of reducing predictive power. Using LASSO's variable screening property (with high probability, the Lasso estimated model includes the substantial covariates and drops the redundant ones) this function can identify the relevant genotype and environment variables in the Variable Methylated Region (VMR) dataset. For each VMR, LASSO is ran three times: 1) Including only the genotype variables for the selection step, 2) Including only the environmental variables for the selection step, and 3) Including both the genotype and environmental variables in the selection step. This is done to ensure that the function captures the variables that are relevant with the variables of their own category (e.g. SNPs that are strongly associated with the DNAme levels of a VMR in the presence of the rest of the SNPs) or in the presence of the covariates of the other category (e.g. SNPs that are strongly associated with the DNAme levels of a VMR in the presence of the rest of BOTH the SNPs AND environmental variables).
#'
#' Each LASSO model uses a tuned lambda that minimizes the 5-fold cross-validation error within its corresponding data. This function uses the lambda.min value in contrast to lambda.1se because its goal within the RAMEN package is to use LASSO to reduce the number of variables that are going to be used next for fitting pairwise interaction models in *lmGE()*. Since at this step variables are being selected base only on main effects, it is preferable to cast a "wider net" and select a slightly higher number of variables that could potentially have a strong interaction effect when paired with another variable. Furthermore, since LASSO is being used in this case as a screening procedure to select variables that will be fit separately in independent models and compared, the overfitting issue of using lambda.min does not impose a big concern. After finding the best lambda value, the sequence of models is fit by coordinate descent using *glmnet()*.
#'
#' @param VMR_df A data frame converted from a Genomic Ranges object. Recommended to use the output of *RAMEN::findCisSNPs()*. Must have one VMR per row, and contain the following columns: VMR_index" (an unique ID for each VMR in VMR_df AS CHARACTERS) and "SNP" (a column with a list as observation, containing the name of the SNPs surrounding the corresponding VMR).  The SNPs contained in the "SNP" column must be present in the object that is indicated in the genotype_matrix argument, and it must contain all the VMRs contained in summarized_methyl_VMR. VMRs with no surrounding SNPs must have an empty list in the SNP column (either list(NULL), list(NA) or list("")).
#' @param environmental_matrix A matrix of environmental variables. Only numeric values are supported. In case of having factors, it is recommended to encode them as numbers or re-code them into dummy variables if there are more than 2 levels.Each column should correspond to an environmental variable and each row to an individual. The ROW names must be the individual IDs.
#' @param genotype_matrix A matrix of number-encoded genotypes. Each column should correspond to a sample, and each row to a SNP site. The suggested encoding is 2 (AA), 1 (AB) and 0 (BB). The COLUMN names must correspond with individual IDs.
#' @param covariates A matrix containing the variables that will be adjusted for in the final GxE models (e.g. cell type proportions, age, etc.)-.Each column should correspond to a covariate and each row to an individual. The ROW names must be the individual IDs.
#' @param summarized_methyl_VMR A data frame containing each individual's VMR summarized region methylation. It is suggested to use the output of RAMEN::summarizeVMRs(). It must have individuals as rows, and VMRs as columns. The names of the columns must correspond to the index of said VMR, and it must match the index of VMR_df$VMR_index. The names of the ROWS must be the sample IDs, and must match with the IDs of the other matrices.
#' @param seed An integer number that initializes a pseudo-random number generator. Random numbers in this function are created during the lambda cross validation and the LASSO stages.
#'
#' @return A data frame with three columns: VMR_index, selected_genot and selected_env. The two last columns contain lists with a vector inside containing the variables that were selected using the LASSO strategy.
#'
#' @importFrom doRNG %dorng%
#' @export


selectVariables = function(VMR_df,
                           genotype_matrix,
                           environmental_matrix,
                           covariates = NULL,
                           summarized_methyl_VMR,
                           seed = NULL) {
  ## Arguments check
  # Check that genotype_matrix, environmental_matrix, covariate matrix (in case it is provided) and summarized_methyl_VMR have the same samples
  if(!all(rownames(summarized_methyl_VMR) %in% colnames(genotype_matrix))){stop("Individual IDs in summarized_methyl_VMR do not match individual IDs in genotype_matrix")}
  if (!all(rownames(summarized_methyl_VMR) %in% rownames(environmental_matrix))) {stop("Individual IDs in summarized_methyl_VMR do not match individual IDs in environmental_matrix")}
  if(!is.null(covariates)){
    if (!all(rownames(summarized_methyl_VMR) %in% rownames(covariates))){stop("Individual IDs in summarized_methyl_VMR do not match individual IDs in the covariates matrix")}}
  #Check that VMR_df has index and SNP column
  if(!all(c("VMR_index","SNP") %in% colnames(VMR_df))){stop("Please make sure the VMR data frame (VMR_df) contains the columns 'SNP' and 'VMR_index'.")}
  #Check that the SNP column on VMRs_df is a list
  if(!is.list(VMR_df$SNP)){stop("Please make sure the 'SNP' column in VMRs_df is a column of lists")}
  if(!is.character(VMR_df$VMR_index)){stop("Please make sure the 'VMR_index' column in VMRs_df is a column of characters")}
  #Check that genotype, environment and covariates are matrices
  if (!is.matrix(genotype_matrix)){stop("Please make sure the genotype data is provided as a matrix.")}
  if (!is.matrix(environmental_matrix)){stop("Please make sure the environmental data is provided as a matrix.")}
  if (!is.null(covariates)){
    if (!is.matrix(covariates)){ stop("Please make sure the covariates data is provided as a matrix.")}}

  ## Set the seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  lasso_results = foreach::foreach(VMR_i = iterators::iter(VMR_df, by = "row"), .combine = "rbind") %dorng%{
    #Select summarized VMR information
    summVMRi = summarized_methyl_VMR %>%
      dplyr::select(VMR_i$VMR_index)
    ## Prepare data
    #subset the genotyping data and match genotype, environment and DNAme IDs
    if(VMR_i$SNP %in% list(NULL) | # Catch VMRs with no surrounding SNPs
       VMR_i$SNP %in% list("") |
       VMR_i$SNP %in% list(NA) ){
      genot_VMRi = c()
    } else {
      genot_VMRi = genotype_matrix[unlist(VMR_i$SNP), rownames(summVMRi)] %>%
        t()
    }
    environ_VMRi = environmental_matrix[rownames(summVMRi),]
    environ_genot_VMRi = cbind(genot_VMRi, environ_VMRi)
    #Bind covariates data
    .GlobalEnv$covariates = covariates
    if (!is.null(covariates)){
      covariates_VMRi = covariates[rownames(summVMRi),]
      genot_VMRi = cbind(genot_VMRi,covariates_VMRi)
      environ_VMRi = cbind(environ_VMRi, covariates_VMRi)
      environ_genot_VMRi = cbind(environ_genot_VMRi, covariates_VMRi)
    }

    ### Run LASSOs
    ## Genotype only
    #conduct k-fold cross-validation to find optimal lambda value
    lambda_genot <- glmnet::cv.glmnet(x = genot_VMRi, #Variables
                              y = summVMRi[,VMR_i$VMR_index], #Response
                              alpha = 1,
                              nfolds = 5)$lambda.min
    #Select the variables with a coefficient > 0
    lasso_genot <- glmnet::glmnet(x = genot_VMRi, #Variables
                          y = summVMRi[,VMR_i$VMR_index], #Response
                          alpha = 1,
                          lambda = lambda_genot)
    coef_genot = stats::coef(lasso_genot)
    coef_genot = coef_genot[coef_genot[,1] > 0,]
    selected_vars_genot = names(coef_genot)[-1]
    selected_vars_genot = selected_vars_genot[!selected_vars_genot %in% colnames(covariates)] #Remove covariates in case they are selected

    #Environment only
    #conduct k-fold cross-validation to find optimal lambda value
    lambda_env <- glmnet::cv.glmnet(x = environ_VMRi, #Variables
                            y = summVMRi[,VMR_i$VMR_index], #Response
                            alpha = 1,
                            nfolds = 5)$lambda.min
    #Select the variables with a coefficient > 0
    lasso_env <- glmnet::glmnet(x = environ_VMRi, #Variables
                        y = summVMRi[,VMR_i$VMR_index], #Response
                        alpha = 1,
                        lambda = lambda_env)
    coef_env = stats::coef(lasso_env)
    coef_env = coef_env[coef_env[,1]>0,]
    selected_vars_env = names(coef_env)[-1] #Remove the intercept from the variables
    selected_vars_env = selected_vars_env[!selected_vars_env %in% colnames(covariates)] #Remove covariates in case they are selected

    #Joint (environment + genotype)
    #conduct k-fold cross-validation to find optimal lambda value
    lambda_joint <- glmnet::cv.glmnet(x = environ_genot_VMRi, #Variables
                              y = summVMRi[,VMR_i$VMR_index], #Response
                              alpha = 1,
                              nfolds = 5)$lambda.min
    #Select the variables with a coefficient > 0
    lasso_joint <- glmnet::glmnet(x = environ_genot_VMRi, #Variables
                          y = summVMRi[,VMR_i$VMR_index], #Response
                          alpha = 1,
                          lambda = lambda_joint)
    coef_joint = stats::coef(lasso_joint)
    coef_joint = coef_joint[coef_joint[,1]>0,]
    selected_vars_joint = names(coef_joint)[-1] #Remove the intercept from the variables
    selected_vars_joint = selected_vars_joint[!selected_vars_joint %in% colnames(covariates)] #Remove covariates in case they are selected

    #Merge results
    selected_union_genot = c(selected_vars_genot, selected_vars_joint) %>%
      unique() %>%
      dplyr::setdiff(colnames(environ_VMRi)) #Remove environmental variables and covariates from the joint selection
    selected_union_env = c(selected_vars_env, selected_vars_joint) %>%
      unique() %>%
      dplyr::setdiff(colnames(genot_VMRi)) #Remove genotype variables and covariates from the joint selection

    ### Create final data frame
    selected_variables_final = data.frame(
      VMR_index = VMR_i$VMR_index)
    selected_variables_final$selected_genot = list(selected_union_genot)
    selected_variables_final$selected_env = list(selected_union_env)
    selected_variables_final
  }
  return(lasso_results)
}
