#' Simulate a delta R squared null distribution of G and E effects on DNAme variability
#'
#' This function simulates the delta R squared distribution under the null hypothesis of G and E having no association with DNA methylation (DNAme) variability through a permutation analysis. To do so, this function shuffles the G and E variables in the dataset, which is followed by a the variable selection and modelling steps with *selectVariables()* and *lmGE()*.These steps are repeated several times as indicated in the *permutations* parameter. By using shuffled G and E data, we simulate the increase of R2 that would be observed in random data using the RAMEN methodology.
#'
#' The core pipeline from the RAMEN package identifies the best explanatory model per VMR. However, despite these models being winners in comparison to models including any other G/E variable(s) in the dataset, some winning models might perform no better than what we would expect by chance. Therefore, the goal of this function is to create a distribution of increase in R2 under the null hypothesis of G and E having no associations with DNAme. The null distribution is obtained through shuffling the G and E variables in a given dataset and conducting the variable selection and G/E model selection. That way, we can simulate how much additional variance would be explained by the models defined as winners by the RAMEN methodology in a scenario where the G and E associations with DNAme are randomized. This distribution can be then used to filter out winning models in the non-shuffled dataset that do not add more to the explained variance of the basal model than what randomized data do.
#'
#' Under the assumption that after adjusting for the concomitant variables all VMRs across the genome follow the same behavior regarding an increment of explained variance with randomized G and E data, we can pool the delta R squared values from all VMRs to create a null distribution taking advantage of the high number of VMRs in the dataset. This assumption decreases significantly the number of permutations required to create a null distribution and reduces the computational time. For further information please read the RAMEN paper (in preparation).
#'
#' @param permutations description
#' @inheritParams selectVariables
#' @inheritParams lmGE
#'
#' @return A data frame with the following columns:
#'  - VMR_index: The unique ID of the VMR.
#'  - model_group: The group to which the winning model belongs to (i.e., G, E, G+E or GxE)
#'  - tot_r_squared: R squared of the winning model
#'  - R2_difference: the increase in R squared obtained by including the G/E variable(s) from the winning model (i.e., the R squared difference between the winning model and the model only with the concomitant variables specified in *covariates*; tot_r_squared - basal_rsquared in the lmGE output)
#'  - AIC_difference/BIC_difference: the AIC/BIC difference between the winning model and the model only with the concomitant variables specified in *covariates*; BIC/AIC - basal_BIC/basal_BIC in the lmGE output)
#'
#' @importFrom foreach %do%
#' @export
#'

nullDistGE = function(VMRs_df,
                     genotype_matrix,
                     environmental_matrix,
                     summarized_methyl_VMR,
                     permutations = 10,
                     covariates = NULL,
                     seed = NULL,
                     model_selection = "AIC"
){
  #Get the shuffle order
  if (!is.null(seed)) set.seed(seed)
  permutation_order = data.frame(sample(rownames(summarized_methyl_VMR),
                                        size = length(rownames(summarized_methyl_VMR))))
  for (i in 1:(permutations-1)){
    permutation_order= cbind(permutation_order,
                             data.frame(sample(rownames(summarized_methyl_VMR),
                                               size = length(rownames(summarized_methyl_VMR)))))
  }
  colnames(permutation_order) = 1:permutations

  #Put the environmental and genotype matrix in the same order to the summarized VMR object
  genotype_matrix = genotype_matrix[,rownames(summarized_methyl_VMR)]
  environmental_matrix = environmental_matrix[rownames(summarized_methyl_VMR),]

  # Permutation analysis
  null_dist = foreach::foreach(i = 1:permutations, .combine = rbind) %do% {
    #Shuffle the datasets
    permutated_genotype = genotype_matrix[,permutation_order[,i]] %>%
      as.matrix()
    rownames(permutated_genotype) = rownames(genotype_matrix)
    colnames(permutated_genotype) = colnames(genotype_matrix)
    permutated_environment = environmental_matrix[permutation_order[,i],] %>%
      as.matrix()
    colnames(permutated_environment) = colnames(environmental_matrix)
    rownames(permutated_environment) = rownames(environmental_matrix)

    # Run RAMEN
    selected_variables = RAMEN::selectVariables(VMRs_df = VMRs_df,
                                                genotype_matrix = permutated_genotype,
                                                environmental_matrix = permutated_environment,
                                                covariates = covariates,
                                                summarized_methyl_VMR = summarized_methyl_VMR,
                                                seed = 1)

    lmGE_res = RAMEN::lmGE(selected_variables = selected_variables,
                           summarized_methyl_VMR = summarized_methyl_VMR,
                           genotype_matrix = permutated_genotype,
                           environmental_matrix = permutated_environment,
                           covariates = covariates,
                           model_selection = model_selection)
    if(model_selection=="AIC"){
      results_perm = data.frame(VMR_index = lmGE_res$VMR_index,
                                tot_r_squared = lmGE_res$tot_r_squared,
                                model_group = lmGE_res$model_group,
                                R2_difference = lmGE_res$tot_r_squared - lmGE_res$basal_rsquared,
                                AIC_difference = lmGE_res$AIC - lmGE_res$basal_rsquared)
    } else if (model_selection=="BIC"){
      results_perm = data.frame(VMR_index = lmGE_res$VMR_index,
                                model_group = lmGE_res$model_group,
                                tot_r_squared = lmGE_res$tot_r_squared,
                                R2_difference = lmGE_res$tot_r_squared - lmGE_res$basal_rsquared,
                                BIC_difference = lmGE_res$BIC - lmGE_res$basal_rsquared)
    }
    results_perm$permutation = i #add the number of permutation
    results_perm
  }

  return(null_dist)
}

