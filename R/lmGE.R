#' Compute G/E/G+E/GxE linear models for a VMRs
#'
#' This will be coupled with apply(VMR_df, 1, FUN = g_e_modelling,...)
#' metadata and environmental objects must have a column named "SubjectNumber" that match the individuals
#'
#' @param VMR ?a single VMR that comes from a VMR_df object. This VMR_df must have ?
#' @param methylation_data ?A data frame where each row is a probe and each column is an individual.
#' The names of the individuals (i.e. colnames) must be the same as the "SubjectNumber" in the metadata
#' and environmental objects
#' @param environmental_data A data frame with one sample per row. The first column must be the "SubjectNumber",
#' and the second one a single environmental variable of choice.
#' @param genotyping_data A data frame with one sample per row. It must contain the columns "SubjectNumber"
#' and "genotype_encoded", which corresponds to the genotypes recoded into a count of 0, 1 or 2 (as characters)
#' representing the number of minor allele copies.
#' @param covariates A data frame with all the covariate information. It must have a "SubjectNumber" column with
#' the individual ID, which matches with the column in genotyping data and environmental data.
#' @param formula_covariates A string containing the covariates to be used in the models separated by a "+"
#' sign (e.g., Sex + Sample_Age_Days + counts.CD4T + counts.NK + counts.Bcell + ...)
#'
#' @return ?a list with the following elements:
#'  - VMR_index
#'  - tagCpG
#'  - SNP
#'  - environment
#'  - statistics
#'  - winning_model
#' @export

lmGE = function(VMR, #VMRs_df_with_cisSNPs
                methylation_data,#methylation data
                environmental_data,
                genotyping_data,
                covariates,
                formula_covariates){

  #Step 1: create a dataframe that contains all the information
  condensated_df =
    #### Query the methylation data
    methylation_data %>%
    tibble::rownames_to_column(var = "probe") %>%
    dplyr::filter(probe %in% VMR["probes"]) %>%
    dplyr::mutate(probe = "tagCpG") %>%
    tibble::column_to_rownames(var = "probe") %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "SubjectNumber") %>%
    #### Query the genotyping data
    dplyr:: left_join(genotyping_data %>%
                        dplyr::filter(SNP.Name == VMR["SNP"]),
                      dplyr::select(SubjectNumber, genotype_encoded),
                      by = "SubjectNumber") %>%
    #### Query the environmental data
    dplyr::left_join(environmental_data,
                     by = "SubjectNumber") %>%
    #### Query the covariates
    dplyr::left_join(covariates,
                     by = "SubjectNumber")
  #Model G
  g_lm = stats::lm(data = condensated_df, formula = stringr::str_glue("tagCpG ~ genotype_encoded + ", formula_covariates) )

  #Model E
  e_lm = stats::lm(data = condensated_df, formula = stringr::str_glue("tagCpG ~ ", colnames(environmental_data[2]), " + ", formula_covariates) )

  #Model G + E
  g_e_lm = stats::lm(data = condensated_df, formula = stringr::str_glue("tagCpG ~ genotype_encoded + ", colnames(environmental_data[2]), " + ", formula_covariates) )

  #Model G x E
  gxe_lm = stats::lm(data = condensated_df, formula = stringr::str_glue("tagCpG ~ genotype_encoded + ", colnames(environmental_data[2]), " +  genotype_encoded*",colnames(environmental_data[2]), " + ", formula_covariates) )

  #AIC and other statistics
  results = stats::AIC(g_lm, e_lm, g_e_lm, gxe_lm) %>%
    dplyr::mutate(r.squared = c(summary(g_lm)$r.squared,
                                summary(e_lm)$r.squared,
                                summary(g_e_lm)$r.squared,
                                summary(gxe_lm)$r.squared),
                  adj.r.squared = c(summary(g_lm)$adj.r.squared,
                                    summary(e_lm)$adj.r.squared,
                                    summary(g_e_lm)$adj.r.squared,
                                    summary(gxe_lm)$adj.r.squared))
  winning_model = rownames(results %>% dplyr::arrange(AIC))[1]

  return(list(VMR_index = VMR["VMR_index"],
              tagCpG = VMR["probes"],
              SNP = VMR["SNP"],
              environment = colnames(environmental_data)[2],
              statistics = results,
              winning_model = winning_model
              #Pvalue of the F test comparing the tests
  ))
}
