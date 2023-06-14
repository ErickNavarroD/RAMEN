#' Summarize Varaible Methylated Rregion's methylation state
#'
#' For each VMR in a dataset, returns an object with the representative methylation score of that region per individual.
#' OPTIONAL: This function supports parallel computing for increased speed. To do so, you have to set the parallel backend
#' in your R session BEFORE running the function (e.g., doFuture::registerDoFuture()) and then the evaluation strategy (e.g., future::plan(multisession)). After that,
#' the function can be run as usual.
#'
#' @param VMRs_df A GRanges object converted to a data frame. Must contain the following columns:
#' "seqnames", "start", "end"  (all of which are produced automatically when doing the object conversion)
#' and "probes" (containing a list where each element contains a vector with the probes
#' constituting the VMR)
#' @param methylation_data A data frame containing M or B values, with samples as columns and probes as rows
#'
#' @return A data frame with samples as rows, and VMRs as columns. The value inside each cell corresponds to the summarized
#' methylation value of said VMR in said individual. The column names correspond to the VMR_index, which is created if not
#' already existing based on the rownames of the VMR_df.
#'
#' @importFrom foreach %dopar%
#' @export

summarizeVMRs = function(VMRs_df,
                         methylation_data){
  if(!"VMR_index" %in% colnames(VMRs_df)){ # Add a VMR index to each region if not already existing
    VMRs_df = VMRs_df %>%
      tibble::rownames_to_column(var = "VMR_index")
  }

  # Check that probes is a list.
  if(class(VMRs_df$probes) != "list"){
    stop("Please make sure the 'probes' column in VMRs_df is a column of lists")
  }

  summarized_VMRs = foreach::foreach(i = VMRs_df$VMR_index, .combine = "cbind") %dopar% {
    probes = VMRs_df %>%
      dplyr::filter(VMR_index == i) %>%
      dplyr::pull(probes) %>%
      unlist()
    subset_meth =  methylation_data[probes, ] %>%
      t() %>%
      as.data.frame()
    median = data.frame(apply(subset_meth,1,median))
    colnames(median) = i
    median
  }
  return(summarized_VMRs)
}
