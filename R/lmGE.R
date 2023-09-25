
#' Fit linear G, E, G+E and GxE models and select the winning model
#'
#' For a set of Variable Methylated Region (VMR), this function fits a set of genotype (G), environment (E), pairwise additive (G + E) or pairwise interaction (G x E) models, one variable at a time, and selects the best fitting one. Additional information for each winning model is provided, such as its R2, its significance from an F-test comparing it to a basal model (i.e., only the covariates), the delta AIC/BIC to the next best model from a different category, and the explained variance decomposed for the G, E and GxE components (when applicable).
#'
#' This function supports parallel computing for increased speed. To do so, you have to set the parallel backend
#' in your R session before running the function (e.g., doFuture::registerDoFuture()) and then the evaluation strategy (e.g., future::plan(multisession)). After that,
#' the function can be run as usual. It is recommended to also set options(future.globals.maxSize= +Inf).
#'
#' For each VMR, this function computes a set of models using the variables indicated in the selected_variables object. From the indicated G and E variables, lmGE() fits four groups of models:
#'  - G: Genetics model - fitted one SNP at a time.
#'  - E: Environmental model - fitted one environmental variable at a time.
#'  - G+E: Additive model - fitted for each pairwise combination of G and E variables indicated in selected_variables.
#'  - GxE: Interaction model - fitted for each pairwise combination of G and E variables indicated in selected_variables.
#'
#'  These models are fit only if the VMR has G or E variables in the selected_variables object. If a VMR does not have neither G nor E variables, that VMR will be ignored and will not be returned in the output object.
#'
#' **Model selection**
#'
#' Following the model fitting stage, the best model **per group** is selected using Akaike Information Criterion (AIC) or Bayesian Information Criterion (BIC). Both of these metrics are statistical approaches to select the best model in the same data set. Both of these metrics have strengths and limitations that make them excel in different situations. We recommend using AIC because BIC assumes that the true model is in the set of compared models. Since this function fits models with individual variables, and we assume that DNAme variability is more likely to be influenced by more than one single SNP/environmental exposure at a time, we hypothesize that in most cases, the true model will not be in the set of compared models. Also, AIC excels in situations where all models in the model space are "incorrect", and AIC is preferentially used in cases where the true underlying function is unknown and our selected model could belong to a very large class of functions where the relationship could be pretty complex. It is worth mentioning however that, both metrics tend to pick the same model in a large number of scenarios. We suggest the users to read Arijit Chakrabarti & Jayanta K. Ghosh, 2011 for further information about the difference between these metrics.
#'
#' After selecting the best model per group (G,E,G+E pr GxE), the model with the lowest AIC or BIC is declared as the winning model. The delta AIC/BIC and difference of R2 is computed relative to the model with the second lowest AIC/BIC (i.e., the best model from a different group to the winning one), and reported in the final object.
#'
#' **Analysis of variance and variance decomposition**
#'
#' After defining the best explanatory model for each VMR, an ANOVA test is conducted to test whether this model is significantly better at explaining the DNAme variability of the respective region than the basal model (i.e., the model including only the specified covariates). The respective F- and p-values are reported in the final returned object. Additionally, a False Discovery Rate (FDR) method is applied to the p-values to address the multiple hypothesis testing.
#' Finally, the variance is decomposed and the relative R2 contribution of each of the variables of interest (G, E and GxE) is reported. This decomposition is done using the relaimpo R package, using the Lindeman, Merenda and Gold (lmg) method, which is based on the heuristic approach of averaging the relative R contribution of each variable over all input orders in the linear model. The estimation of the partitioned R2 of each factor in the models was conducted keeping the covariates always in the model as first entry (i.e., the variables specified in covariates did not change order). For further information, we suggest the users to read the documentation and publication of the relaimpo R package (Grömping, 2006).
#'
#' @param selected_variables A data frame obtained with *RAMEN::selectVariables()*. This data frame must contain three columns: 'VMR_index' with characters of an unique ID of each VMR; ´selected_genot' and 'selected_env' with the SNPs and environmental variables, respectively, that will be used for fitting the genotype (G), environment (E), additive (G + E) or interaction (G x E) models. The columns 'selected_env' and 'selected_genot' must contain lists as elements; VMRs with no environmental or genotype selected variables must contain an empty list with NULL, NA , character(0) or "" inside.
#' @param environmental_matrix A matrix of environmental variables. Only numeric values are supported. In case of factor variables, it is recommended to encode them as numbers or re-code them into dummy variables if there are more than two levels. Columns must correspond to environmental variables and rows to individuals. Row names must be the individual IDs.
#' @param genotype_matrix A matrix of number-encoded genotypes. Columns must correspond to samples, and rows to SNPs. We suggest using a gene-dosage model, which would encode the SNPs ordinally depending on the genotype allele charge, such as 2 (AA), 1 (AB) and 0 (BB). The column names must correspond with individual IDs.
#' @param summarized_methyl_VMR A data frame containing each individual's VMR summarized region methylation. It is suggested to use the output of RAMEN::summarizeVMRs().Rows must reflects individuals, and columns VMRs The names of the columns must correspond to the index of said VMR, and it must match the index of VMR_df$VMR_index. The names of the rows must correspond to the sample IDs, and must match with the IDs of the other matrices.
#' @param covariates A matrix containing the concomitant variables that will always be adjusted for in the final G/E/G+E/GxE models (e.g. cell type proportions, age, etc.).Column must correspond to covariates and rows to individuals. Row names must be the individual IDs.
#' @param model_selection Which metric to use to select the best model for each VMR. Supported options are "AIC" or BIC". More information about which one to use can be found in the Details section.
#' @return A data frame with the following columns:
#'  - VMR_index: The unique ID of the VMR
#'  - model_group: The group to which the winning model belongs to (i.e., G, E, G+E or GxE)
#'  - variables: The variable(s) that are present in the winning model (excluding the covariates, which are included in all the models)
#'  - tot_r_squared: R squared of the winning model
#'  - g_r_squared: Estimated R2 allocated to the G in the winning model, if applicable.
#'  - e_r_squared: Estimated R2 allocated to the E in the winning model, if applicable.
#'  - gxe_r_squared: Estimated R2 allocated to the interaction in the winning model (GxE), if applicable.
#'  - AIC/BIC: Model comparison value used to select the best model in each VMR (depending on the option specified in the argument model_selection).
#'  - second_winner: The second group that possesses the next best model after the winning one (i.e., G, E, G+E or GxE). This column may have NA if the variables in selected_variables correspond only to one group (G or E), so that there is no other groups to compare to.
#'  - delta_aic/delta_bic: The difference of AIC or BIC value (depending on the option specified in the argument model_selection) of the winning model and the best model from the second_winner group (i.e., G, E, G+E or GxE). This column may have NA if the variables in selected_variables correspond only to one group (G or E), so that there is no other groups to compare to.
#'  - delta_r_squared: The R2 of the winning model - R2 of the second winner model. This column may have NA if the variables in selected_variables correspond only to one group (G or E), so that there is no other groups to compare to.
#'  - basal_AIC/basal_BIC: AIC or BIC of the basal model (i.e., model fitted only with the concomitant variables specified in the *covariates* argument)
#'  - basal_rsquared: The R2 of the basal model (i.e., model fitted only with the concomitant variables specified in the *covariates* argument)
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @export
#'
lmGE = function(selected_variables,
                  summarized_methyl_VMR,
                  genotype_matrix,
                  environmental_matrix,
                  covariates = NULL,
                  model_selection = "AIC"){
  #Check arguments
  # Check that genotype_matrix, environmental_matrix, covariate matrix (in case it is provided) and summarized_methyl_VMR have the same samples
  if(!all(rownames(summarized_methyl_VMR) %in% colnames(genotype_matrix))) stop("Individual IDs in summarized_methyl_VMR do not match individual IDs in genotype_matrix")
  if (!all(rownames(summarized_methyl_VMR) %in% rownames(environmental_matrix))) stop("Individual IDs in summarized_methyl_VMR do not match individual IDs in environmental_matrix")
  if(!is.null(covariates)){
    if (!all(rownames(summarized_methyl_VMR) %in% rownames(covariates)))stop("Individual IDs in summarized_methyl_VMR do not match individual IDs in the covariates matrix")}
  #Check that selected_variables has the right columns
  if(!all(c("VMR_index","selected_genot", "selected_env") %in% colnames(selected_variables))) stop("Please make sure the selected_variables data frame contains the columns 'VMR_index', 'selected_genot' and 'selected_env'.")
  #Check that the selected_genot and selected_env columns on selected_variables is a list and the index is characters
  if(!is.list(selected_variables$selected_genot)) stop("Please make sure the 'selected_genot' column in selected_variables contains lists as elements")
  if(!is.list(selected_variables$selected_env)) stop("Please make sure the 'selected_env' column in selected_variables contains lists as elements")
  if(!is.character(selected_variables$VMR_index)) stop("Please make sure the 'VMR_index' column in selected_variables contains characters")
  #Check that genotype, environment and covariates are matrices
  if (!is.matrix(genotype_matrix)) stop("Please make sure the genotype data is provided as a matrix.")
  if (!is.matrix(environmental_matrix)) stop("Please make sure the environmental data is provided as a matrix.")
  if (!is.null(covariates)){
    if (!is.matrix(covariates)) stop("Please make sure the covariates data is provided as a matrix.")}
  if(!model_selection %in% c("AIC", "BIC")) stop("Please make sure your model_selection method is 'AIC' or 'BIC'")

  #Remove VMRs that have no selected G and no selected E
  selected_variables = selected_variables %>%
    dplyr::filter(!(selected_env %in% c(list(NULL), list(""), list(NA), list(character(0))) &
                     selected_genot %in% c(list(NULL), list(""), list(NA), list(character(0)))))

  #Select the winning model
  winning_models = foreach::foreach(VMR_i = iterators::iter(selected_variables, by = "row"),
                                    .combine = "rbind") %dopar% { #For every VMR
                                      #Create the data frame with all the information for each VMR
                                      summ_vmr_i = as.matrix(summarized_methyl_VMR[,VMR_i$VMR_index])
                                      colnames(summ_vmr_i) = "DNAme"
                                      if (!VMR_i$selected_env %in% c(list(NULL), list(""), list(NA), list(character(0)))) {
                                        if(length(VMR_i$selected_env[[1]]) == 1){
                                          env_i = environmental_matrix[rownames(summarized_methyl_VMR), unlist(VMR_i$selected_env)] %>%
                                            as.matrix()
                                          colnames(env_i) = unlist(VMR_i$selected_env)
                                        } else env_i = environmental_matrix[rownames(summarized_methyl_VMR), unlist(VMR_i$selected_env)]
                                      } else env_i = NULL
                                      if (!VMR_i$selected_genot %in% c(list(NULL), list(""), list(NA), list(character(0)))) {
                                        if(length(VMR_i$selected_genot[[1]]) == 1 ){
                                          genot_i = genotype_matrix[unlist(VMR_i$selected_genot),rownames(summarized_methyl_VMR)] %>%
                                            as.matrix()
                                          colnames(genot_i) = unlist(VMR_i$selected_genot)
                                        } else {
                                          genot_i = genotype_matrix[unlist(VMR_i$selected_genot),rownames(summarized_methyl_VMR)] %>%
                                            t()
                                        }
                                      } else genot_i = NULL
                                      covariates_i = covariates[rownames(summarized_methyl_VMR),]
                                      full_data_vmr_i = cbind(summ_vmr_i, env_i, genot_i, covariates_i)
                                      colnames(full_data_vmr_i) = make.names(colnames(full_data_vmr_i))
                                      #Set the basal model (only covariates)
                                      basal_model_formula = colnames(covariates) %>%
                                        make.names() %>%
                                        paste( collapse = " + ")

                                      ## Fit models involving G if G has selected variables
                                      if (!VMR_i$selected_genot %in% c(list(NULL), list(""), list(NA), list(character(0)))) {
                                        models_g_involving_df = foreach::foreach(SNP = unlist(VMR_i$selected_genot),
                                                                                 .combine = "rbind") %do% { #For each SNP
                                                                                   ### Fit G models
                                                                                   model_g =  stats::lm(data = as.data.frame(full_data_vmr_i), formula = stringr::str_glue("DNAme ~ ",  make.names(SNP), " + ", basal_model_formula) )

                                                                                   #Create data frame structure for the results
                                                                                   model_g_df = data.frame(model_group = "G")
                                                                                   model_g_df$variables = list(SNP)
                                                                                   if(model_selection == "AIC") model_g_df$AIC = stats::AIC(model_g)
                                                                                   if(model_selection == "BIC") model_g_df$BIC = stats::BIC(model_g)
                                                                                   model_g_df$tot_r_squared = summary(model_g)$r.squared
                                                                                   #model_g_df$tot_adj_r_squared = summary(model_g)$adj.r.squared

                                                                                   if (!VMR_i$selected_env %in% c(list(NULL), list(""), list(NA), list(character(0)))){
                                                                                     ### Fit GxE and G+E models if E is not empty
                                                                                     models_joint_df = foreach::foreach(env = unlist(VMR_i$selected_env),  #For every env var
                                                                                                               .combine = "rbind") %do% {
                                                                                                                 #Fit G + E
                                                                                                                 model_ge = stats::lm(data = as.data.frame(full_data_vmr_i), formula = stringr::str_glue("DNAme ~ ",  make.names(SNP), " + ", make.names(env), " + ", basal_model_formula) )

                                                                                                                 #Create data frame structure for the results
                                                                                                                 model_ge_df = data.frame(model_group = "G+E")
                                                                                                                 model_ge_df$variables = list(c(SNP, env))
                                                                                                                 if(model_selection == "AIC") model_ge_df$AIC = stats::AIC(model_ge)
                                                                                                                 if(model_selection == "BIC") model_ge_df$BIC = stats::BIC(model_ge)
                                                                                                                 model_ge_df$tot_r_squared = summary(model_ge)$r.squared
                                                                                                                 #Fit GxE
                                                                                                                 model_gxe = stats::lm(data = as.data.frame(full_data_vmr_i), formula = stringr::str_glue("DNAme ~ ",  make.names(SNP), " + ", make.names(env), " + ",  make.names(SNP), "*", make.names(env), " + ", basal_model_formula) )

                                                                                                                 #Create data frame structure for the results
                                                                                                                 model_gxe_df = data.frame(model_group = "GxE")
                                                                                                                 model_gxe_df$variables = list(c(SNP, env))
                                                                                                                 if(model_selection == "AIC") model_gxe_df$AIC = stats::AIC(model_gxe)
                                                                                                                 if(model_selection == "BIC") model_gxe_df$BIC = stats::BIC(model_gxe)
                                                                                                                 model_gxe_df$tot_r_squared = summary(model_gxe)$r.squared

                                                                                                                 #Return joint models
                                                                                                                 temp_models_joint = rbind(model_gxe_df, model_ge_df)
                                                                                                                 temp_models_joint
                                                                                                               }
                                                                                   } else models_joint_df = NULL

                                                                                   #Return object with all the G-involved models
                                                                                   temp_models_g_involving = rbind(model_g_df, models_joint_df)
                                                                                   temp_models_g_involving
                                                                                 }
                                      } else models_g_involving_df = NULL

                                      ### Compute E models if E is not empty
                                      if (!VMR_i$selected_env %in% c(list(NULL), list(""), list(NA), list(character(0)))){ #For each env var
                                        models_e_df = foreach::foreach(env = unlist(VMR_i$selected_env),  #For every env var
                                                              .combine = "rbind") %do% {
                                                                #Fit E models
                                                                model_e = stats::lm(data = as.data.frame(full_data_vmr_i), formula = stringr::str_glue("DNAme ~ ", make.names(env), " + ", basal_model_formula) )

                                                                #Create data frame structure for the results
                                                                model_e_df = data.frame(model_group = "E")
                                                                model_e_df$variables = list(c(env))
                                                                if(model_selection == "AIC") model_e_df$AIC = stats::AIC(model_e)
                                                                if(model_selection == "BIC") model_e_df$BIC = stats::BIC(model_e)
                                                                model_e_df$tot_r_squared = summary(model_e)$r.squared
                                                                #Return the final object
                                                                model_e_df
                                                              }
                                      } else models_e_df = NULL

                                      #Create object with the metrics for all the fitted models
                                      all_models_VMR_i = rbind(models_g_involving_df, models_e_df)

                                      #Select the best model per category (G,E,GxE,G+E) and compute its delta AIC/BIC
                                      if(model_selection == "AIC"){
                                        best_models_VMR_i = all_models_VMR_i %>%
                                          dplyr::group_by(model_group) %>%
                                          dplyr::filter(AIC == min(AIC)) %>%
                                          dplyr::slice(1) %>%  #In case there are more than one model per group with the exact same AIC, pick the first one
                                          dplyr::arrange(AIC, dplyr::desc(tot_r_squared)) %>%
                                          dplyr::ungroup() %>%
                                          dplyr::mutate(delta_aic = abs(AIC - dplyr::lead(AIC)))
                                      } else if (model_selection == "BIC"){
                                        best_models_VMR_i = all_models_VMR_i %>%
                                          dplyr::group_by(model_group) %>%
                                          dplyr::filter(BIC == min(BIC)) %>%
                                          dplyr::slice(1) %>% #In case there are more than one model per group with the exact same AIC, pick the first one
                                          dplyr::arrange(BIC,dplyr::desc(tot_r_squared) ) %>%
                                          dplyr::ungroup() %>%
                                          dplyr::mutate(delta_bic = abs(BIC - dplyr::lead(BIC)))
                                      }


                                      #Create the final object that will be returned
                                      if(model_selection == "AIC"){
                                        winning_model_VMR_i = best_models_VMR_i %>%
                                          dplyr::filter(AIC == min(AIC)) %>%
                                          #In case there is more than one model with the exact same AIC from different groups, pick the one with the highest tot_r_squared
                                          dplyr::slice(1) %>%
                                          dplyr::mutate(second_winner = best_models_VMR_i$model_group[2],
                                                        delta_r_squared = best_models_VMR_i$tot_r_squared[1] - best_models_VMR_i$tot_r_squared[2])
                                      }else if(model_selection == "BIC"){
                                        winning_model_VMR_i = best_models_VMR_i %>%
                                          dplyr::filter(BIC == min(BIC)) %>%
                                          #In case there is more than one model with the exact same AIC from different groups, pick the one with the highest tot_r_squared
                                          dplyr::slice(1) %>%
                                          dplyr::mutate(second_winner = best_models_VMR_i$model_group[2],
                                                        delta_r_squared = best_models_VMR_i$tot_r_squared[1] - best_models_VMR_i$tot_r_squared[2])
                                      }

                                      #Test the winning model against the basal one and decompose variance for the G, E and GxE components
                                      model_basal = stats::lm(data = as.data.frame(full_data_vmr_i), formula = stringr::str_glue("DNAme ~ ", basal_model_formula) ) #set the basal model for comparing the rest
                                      if (model_selection == "AIC"){
                                        winning_model_VMR_i$basal_AIC = stats::AIC(model_basal)
                                      } else if(model_selection == "BIC"){
                                        winning_model_VMR_i$basal_BIC = stats::BIC(model_basal)}
                                      winning_model_VMR_i$basal_rsquared = summary(model_basal)$r.squared
                                      if(winning_model_VMR_i$model_group == "G"){
                                        winning_lm = stats::lm(data = as.data.frame(full_data_vmr_i), formula = stringr::str_glue("DNAme ~ ", make.names(unlist(winning_model_VMR_i$variables)), " + ", basal_model_formula) )
                                        r_decomp = relaimpo::calc.relimp.lm(object =  winning_lm ,
                                                                            rela = FALSE,
                                                                            type = "last") #This would be the equivalent to using lmg and setting always = covariates.
                                        winning_model_VMR_i$g_r_squared = r_decomp$last[make.names(unlist(winning_model_VMR_i$variables))[1]]
                                        winning_model_VMR_i$e_r_squared = NA_real_
                                        winning_model_VMR_i$gxe_r_squared = NA_real_
                                      }else if (winning_model_VMR_i$model_group == "E"){
                                        winning_lm = stats::lm(data = as.data.frame(full_data_vmr_i), formula = stringr::str_glue("DNAme ~ ", make.names(unlist(winning_model_VMR_i$variables))[1], " + ", basal_model_formula) )
                                        r_decomp = relaimpo::calc.relimp.lm(object =  winning_lm,
                                                                            rela = FALSE,
                                                                            type = "last") #This would be the equivalent to using lmg and setting always = covariates.
                                        winning_model_VMR_i$g_r_squared = NA_real_
                                        winning_model_VMR_i$e_r_squared = r_decomp$last[make.names(unlist(winning_model_VMR_i$variables))[1]]
                                        winning_model_VMR_i$gxe_r_squared = NA_real_
                                      }else if (winning_model_VMR_i$model_group == "G+E"){
                                        winning_lm = stats::lm(data = as.data.frame(full_data_vmr_i), formula = stringr::str_glue("DNAme ~ ", make.names(unlist(winning_model_VMR_i$variables))[1], " + ",make.names(unlist(winning_model_VMR_i$variables))[2], " + ", basal_model_formula) )
                                        r_decomp = relaimpo::calc.relimp.lm(object =  winning_lm,
                                                                            rela = FALSE,
                                                                            type = "lmg",
                                                                            always = colnames(covariates_i))
                                        winning_model_VMR_i$g_r_squared = r_decomp$lmg[make.names(unlist(winning_model_VMR_i$variables))[1]]
                                        winning_model_VMR_i$e_r_squared = r_decomp$lmg[make.names(unlist(winning_model_VMR_i$variables))[2]]
                                        winning_model_VMR_i$gxe_r_squared = NA_real_
                                      }else if (winning_model_VMR_i$model_group == "GxE"){
                                        winning_lm = stats::lm(data = as.data.frame(full_data_vmr_i), formula = stringr::str_glue("DNAme ~ ", make.names(unlist(winning_model_VMR_i$variables))[1], " + ",make.names(unlist(winning_model_VMR_i$variables))[2], " + ", make.names(unlist(winning_model_VMR_i$variables))[1], "*",make.names(unlist(winning_model_VMR_i$variables))[2], " + ", basal_model_formula) )
                                        r_decomp = relaimpo::calc.relimp.lm(object =  winning_lm,
                                                                            rela = FALSE,
                                                                            type = "lmg",
                                                                            always = colnames(covariates_i)) #This slightly underestimates the relative importance compared to not using the covariates as the basal model, but in the interaction option the computational time is greatly increased if the relative contribution of all the other covariates is also estimated (which we dont look at anyways). So, because of the high dimensional nature of this package, this option will be used.
                                        winning_model_VMR_i$g_r_squared = r_decomp$lmg[make.names(unlist(winning_model_VMR_i$variables))[1]]
                                        winning_model_VMR_i$e_r_squared = r_decomp$lmg[make.names(unlist(winning_model_VMR_i$variables))[2]]
                                        winning_model_VMR_i$gxe_r_squared = r_decomp$lmg[stringr::str_glue(make.names(unlist(winning_model_VMR_i$variables))[1], ":",make.names(unlist(winning_model_VMR_i$variables))[2])]
                                      }

                                      winning_model_VMR_i$VMR_index = VMR_i$VMR_index
                                      #Return final object
                                      winning_model_VMR_i
                                    }

  #Compute FDR and rearrange columns
  if (model_selection == "AIC"){
    winning_models = winning_models %>%
      dplyr::select(VMR_index, model_group, variables, tot_r_squared, g_r_squared, e_r_squared, gxe_r_squared, AIC, second_winner, delta_aic, delta_r_squared, basal_AIC, basal_rsquared) %>%
      as.data.frame()
  } else if (model_selection == "BIC"){
    winning_models = winning_models %>%
      dplyr::select(VMR_index, model_group, variables, tot_r_squared, g_r_squared, e_r_squared, gxe_r_squared, BIC, second_winner, delta_bic, delta_r_squared, basal_BIC, basal_rsquared) %>%
      as.data.frame()
  }

  return(winning_models)
}

