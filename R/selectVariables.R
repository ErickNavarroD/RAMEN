#' Selection of environment and genotype variables for Variable Methylated Regions (VMRs)
#'
#' For each VMR, this function selects genotype and environmental variables using LASSO.
#'
#' This function supports parallel computing for increased speed. To do so, you have to set the parallel back-end
#' in your R session before running the function (e.g., doFuture::registerDoFuture()) and then the evaluation strategy (e.g., future::plan(multisession)). After that, the function can be run as usual. It is recommended to also set options(future.globals.maxSize= +Inf).
#'
#' selectVariables() uses LASSO, which is an embedded variable selection method that penalizes models that are more complex (i.e., that contain more variables) in favor of simpler models (i.e. that contain less variables), but not at the expense of reducing predictive power. Using LASSO's variable screening property (with high probability, the LASSO estimated model includes the substantial covariates and drops the redundant ones) this function selects genotype and environment variables with potential relevance in the Variable Methylated Region (VMR) dataset (see also Bühlmann and van de Geer, 2011). For each VMR, LASSO is run three times: 1) including only the genotype variables for the selection step, 2) including only the environmental variables for the selection step, and 3) Including both the genotype and environmental variables in the selection step. This is done to ensure that the function captures the variables that are relevant within their own category (e.g., SNPs that are strongly associated with the DNAme levels of a VMR in the presence of the rest of the SNPs) or in the presence of the variables of the other category (e.g. SNPs that are strongly associated with the DNAme levels of a VMR in the presence of the rest of BOTH the SNPs AND environmental variables). Every time LASSO is run, the basal covariates (i.e., concomitant variables )indicated in the argument *covariates* are not penalized (i.e., those variables are always included in the models and their coefficients are not subjected to shrinkage). That way, only the most promising E and G variables in the presence of the concomitant variables will be selected.
#'
#' Each LASSO model uses a tuned lambda that minimizes the 5-fold cross-validation error within its corresponding data. This function uses the lambda.min value in contrast to lambda.1se because its goal within the RAMEN package is to use LASSO to reduce the number of variables that are going to be used next for fitting pairwise interaction models in *lmGE()*. Since at this step variables are being selected based only on main effects, it is preferable to cast a "wider net" and select a slightly higher number of variables that could potentially have a strong interaction effect when paired with another variable. Furthermore, since in this case LASSO is being used as a screening procedure to select variables that will be fit separately in independent models and compared, the overfitting issue of using lambda.min does not impose a big concern. After finding the best lambda value, the sequence of models is fit by coordinate descent using *glmnet()*. Random numbers in this function are created during the lambda cross validation and the LASSO stages. Setting a seed is highly encouraged for result reproducibility using the *seed* argument. Please note that setting a seed inside of this function modifies the seed globally (which is R's default behavior).
#'
#' Note: If you want to conduct the variable selection step only in one data set (i.e., only in the genotype), you can set the argument *environmental_matrix = NULL*.
#'
#'
#' @param VMRs_df A data frame converted from a GRanges object. Recommended to use the output of *RAMEN::findCisSNPs()*. Must have one VMR per row, and contain the following columns: "VMR_index" (a unique ID for each VMR in VMRs_df AS CHARACTERS) and "SNP" (a column with a list as observation, containing the name of the SNPs surrounding the corresponding VMR).  The SNPs contained in the "SNP" column must be present in the object that is indicated in the genotype_matrix argument, and it must contain all the VMRs contained in summarized_methyl_VMR. VMRs with no surrounding SNPs must have an empty list in the SNP column (either list(NULL), list(NA), list("") or list(character(0)) ).
#' @param environmental_matrix A matrix of environmental variables. Only numeric values are supported. In case of factor variables, it is recommended to encode them as numbers or re-code them into dummy variables if there are more than two levels. Columns must correspond to environmental variables and rows to individuals. Row names must be the individual IDs.
#' @param genotype_matrix A matrix of number-encoded genotypes. Columns must correspond to samples, and rows to SNPs. We suggest using a gene-dosage model, which would encode the SNPs ordinally depending on the genotype allele charge, such as 2 (AA), 1 (AB) and 0 (BB). The column names must correspond with individual IDs.
#' @param summarized_methyl_VMR A data frame containing each individual's VMR summarized region methylation. It is suggested to use the output of RAMEN::summarizeVMRs().Rows must reflects individuals, and columns VMRs The names of the columns must correspond to the index of said VMR, and it must match the index of VMRs_df$VMR_index. The names of the rows must correspond to the sample IDs, and must match with the IDs of the other matrices.
#' @param covariates A matrix containing the covariates (i.e., concomitant variables / variables that are not the ones you are interested in) that will be adjusted for in the final GxE models (e.g., cell type proportions, age, etc.). Each column should correspond to a covariate and each row to an individual. Row names must correspond to the individual IDs.
#' @param seed An integer number that initializes a pseudo-random number generator. Random numbers in this function are created during the lambda cross validation and the LASSO stages. Setting a seed is highly encouraged for result reproducibility. **Please note that setting a seed in this function modifies the seed globally**.
#'
#' @return A data frame with three columns:
#'  - VMR_index: Unique VMR ID.
#'  - selected_genot: List-containing column with the selected SNPs.
#'  - selected_env: List-containing column with the selected environmental variables.
#'
#' @importFrom doRNG %dorng%
#' @export
#'
selectVariables = function(VMRs_df,
                           genotype_matrix,
                           environmental_matrix,
                           covariates = NULL,
                           summarized_methyl_VMR,
                           seed = NULL) {
  ## Arguments check
  # Check that genotype_matrix, environmental_matrix, covariate matrix (in case it is provided) and summarized_methyl_VMR have the same samples
  if(!all(rownames(summarized_methyl_VMR) %in% colnames(genotype_matrix))) stop("Individual IDs in summarized_methyl_VMR do not match individual IDs in genotype_matrix")
  if (!all(rownames(summarized_methyl_VMR) %in% rownames(environmental_matrix))) stop("Individual IDs in summarized_methyl_VMR do not match individual IDs in environmental_matrix")
  if(!is.null(covariates)){
    if (!all(rownames(summarized_methyl_VMR) %in% rownames(covariates)))stop("Individual IDs in summarized_methyl_VMR do not match individual IDs in the covariates matrix")}
  #Check that VMRs_df has index and SNP column
  if(!all(c("VMR_index","SNP") %in% colnames(VMRs_df))) stop("Please make sure the VMR data frame (VMRs_df) contains the columns 'SNP' and 'VMR_index'.")
  #Check that the SNP column on VMRs_df is a list
  if(!is.list(VMRs_df$SNP)) stop("Please make sure the 'SNP' column in VMRs_df is a column of lists")
  if(!is.character(VMRs_df$VMR_index)) stop("Please make sure the 'VMR_index' column in VMRs_df is a column of characters")
  #Check that genotype, environment and covariates are matrices
  if (!is.matrix(genotype_matrix)) stop("Please make sure the genotype data is provided as a matrix.")
  if(!is.null(environmental_matrix)){
    if (!is.matrix(environmental_matrix)) stop("Please make sure the environmental data is provided as a matrix.")}
  if (!is.null(covariates)){
    if (!is.matrix(covariates)) stop("Please make sure the covariates data is provided as a matrix.")}

  ## Set the seed
  if (!is.null(seed)) set.seed(seed)

  lasso_results = foreach::foreach(VMR_i = iterators::iter(VMRs_df, by = "row"), .combine = "rbind") %dorng%{
    #Select summarized VMR information
    summVMRi = summarized_methyl_VMR %>%
      dplyr::select(VMR_i$VMR_index)
    ## Prepare data
    #subset the genotyping data and match genotype, environment and DNAme IDs
    if(VMR_i$SNP %in% list(NULL) | # Catch VMRs with no surrounding SNPs
       VMR_i$SNP %in% list("") |
       VMR_i$SNP %in% list(NA) |
       VMR_i$SNP %in% list(character(0))){
      genot_VMRi = c()
      any_snp = FALSE
    } else if (length(VMR_i$SNP[[1]]) == 1){ #Special case of sub-setting if SNP is only one because the result is a vector and not a matrix
      genot_VMRi = genotype_matrix[unlist(VMR_i$SNP), rownames(summVMRi)] %>%
        as.matrix()
      colnames(genot_VMRi) = VMR_i$SNP[[1]]
      any_snp = TRUE
    } else {
      genot_VMRi = genotype_matrix[unlist(VMR_i$SNP), rownames(summVMRi)] %>%
        t()
      any_snp = TRUE
    }
    if(ncol(environmental_matrix) == 1){
      environ_VMRi = environmental_matrix[rownames(summVMRi),] %>%
        as.matrix()
      colnames(environ_VMRi) = colnames(environmental_matrix)
    } else environ_VMRi = environmental_matrix[rownames(summVMRi),]
    environ_genot_VMRi = cbind(genot_VMRi, environ_VMRi)
    #Bind covariates data
    if (!is.null(covariates)){
      if (ncol(covariates) == 1){
        covariates_VMRi = covariates[rownames(summVMRi),] %>% #Match the covariates dataset with the VMRs information
          as.matrix()
        colnames(covariates_VMRi) = colnames(covariates)
      } else covariates_VMRi = covariates[rownames(summVMRi),]
      genot_VMRi = cbind(genot_VMRi,covariates_VMRi)
      environ_VMRi = cbind(environ_VMRi, covariates_VMRi)
      environ_genot_VMRi = cbind(environ_genot_VMRi, covariates_VMRi)
      ncol_covariates = ncol(covariates_VMRi)
    } else ncol_covariates = 0

    ### Run LASSOs
    ## Genotype only
    #Get coefficients with the optimal lambda found by k-fold cross-validation
    if (any_snp){ #Catch cases when VMRs dont have surrounding genotyped SNPs
      coef_genot = stats::coef(glmnet::cv.glmnet(x = genot_VMRi, #Variables
                                                 y = summVMRi[,VMR_i$VMR_index], #Response
                                                 alpha = 1,
                                                 nfolds = 5,
                                                 penalty.factor = c(rep(1, ncol(genot_VMRi)- ncol_covariates),
                                                                    rep(0, ncol_covariates))), #Unpenalize the variables in covariates (i.e., force LASSO to keep them in all the situations)
                               s = "lambda.min")
      #Select the variables with a coefficient > 0
      coef_genot = coef_genot[abs(coef_genot[,1]) > 0,]
      selected_vars_genot = names(coef_genot)[-1]
      selected_vars_genot = selected_vars_genot[!selected_vars_genot %in% colnames(covariates)] #Remove covariates from selected variables
    } else selected_vars_genot = character(0)

    #Environment only
    #Get coefficients with the optimal lambda found by k-fold cross-validation
    if (!is.null(environ_VMRi)){ #catch scenario where users would not add environmental variables
      coef_env = stats::coef(glmnet::cv.glmnet(x = environ_VMRi, #Variables
                                               y = summVMRi[,VMR_i$VMR_index], #Response
                                               alpha = 1,
                                               nfolds = 5,
                                               penalty.factor = c(rep(1, ncol(environ_VMRi)- ncol_covariates), #Unpenalize the variables in covariates (i.e., force LASSO to keep them in all the situations)
                                                                  rep(0, ncol_covariates))),
                             s = "lambda.min")
      #Select the variables with a coefficient > 0
      coef_env = coef_env[abs(coef_env[,1]) > 0,]
      selected_vars_env = names(coef_env)[-1] #Remove the intercept from the variables
      selected_vars_env = selected_vars_env[!selected_vars_env %in% colnames(covariates)] #Remove covariates from selected variables

      if (any_snp){
        #Joint (environment + genotype) only when we have Genotype and Environmental variables.
        #Get coefficients with the optimal lambda found by k-fold cross-validation
        coef_joint = stats::coef(glmnet::cv.glmnet(x = environ_genot_VMRi, #Variables
                                                   y = summVMRi[,VMR_i$VMR_index], #Response
                                                   alpha = 1,
                                                   nfolds = 5,
                                                   penalty.factor = c(rep(1, ncol(environ_genot_VMRi) - ncol_covariates),
                                                                      rep(0, ncol_covariates))), #Unpenalize the variables in covariates (i.e., force LASSO to keep them in all the situations)
                                 s = "lambda.min")
        #Select the variables with an abs(coefficient) > 0
        coef_joint = coef_joint[abs(coef_joint[,1]) > 0,]
        selected_vars_joint = names(coef_joint)[-1] #Remove the intercept from the variables
        selected_vars_joint = selected_vars_joint[!selected_vars_joint %in% colnames(covariates)] #Remove covariates from selected variables
      } else selected_vars_joint = character(0)
    } else {
      selected_vars_env = character(0)
      selected_vars_joint = character(0)
    }

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
