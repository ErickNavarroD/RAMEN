#' Summarize the methylation states of Variable Methylated Loci (VML)
#'
#' This function computes a representative methylation score for each Variable Methylated Locus (VML) in a dataset. It returns a data frame with the median methylation of each region per individual.
#' For each VML in a dataset, returns a  with the median methylation of that region (columns) per individual (rows) as representative score.
#'
#' This function supports parallel computing for increased speed. To do so, you have to set the parallel backend in your R session BEFORE running the function (e.g., *doParallel::registerDoParallel(4)*). After that,
#' the function can be run as usual.
#'
#' @param VML A GRanges-like data frame. Must contain the following columns:
#' "seqnames", "start", "end" and "probes" (containing lists as elements, where each contains a vector with the probes constituting the VML). This is the "VML" object returned by the *findVML()* function.
#' @param methylation_data A data frame containing M or B values, with samples as columns and probes as rows. Row names must be the CpG probe IDs.
#'
#' @return A data frame with samples as rows, and VML as columns. The value inside each cell corresponds to the summarized methylation value of said VML in the corresponding individual. The column names correspond to the VML_index.
#'
#' @importFrom foreach %dopar%
#' @export

summarizeVML = function(VML,
                         methylation_data){
  if(!"VML_index" %in% colnames(VML)){ # Add a VML index to each region if not already existing
    VML = VML %>%
      tibble::rownames_to_column(var = "VML_index")
  }

  if(!all(unique(unlist(VML$probes)) %in% rownames(methylation_data))){
    warning("Some probes listed in the VML data frame are not found in the methylation data. Please check that all probes listed in the 'probes' column of the VML data frame are present in the row names of the methylation data frame to avoid having NAs.")
  }
  if(!is.data.frame(methylation_data) ){
    if(is.matrix(methylation_data)){
      methylation_data = as.data.frame(methylation_data)
    } else {
      stop("Please make sure the methylation data is a data frame or matrix with samples as columns and probes as rows.")
    }
  }

  # Check that probes is a list.
  if(!is.list(VML$probes)){
    stop("Please make sure the 'probes' column in the VML data frame is a column of lists")
  }

  summarized_VML = foreach::foreach(i = VML$VML_index, .combine = "cbind") %dopar% {
    probes = VML %>%
      dplyr::filter(VML_index == i) %>%
      dplyr::pull(probes) %>%
      unlist()
    subset_meth =  methylation_data[probes, ] %>%
      t() %>%
      as.data.frame()
    median = data.frame(apply(subset_meth,1,median))
    colnames(median) = i
    median
  }
  return(summarized_VML)
}
