#' Compute the median probe methylation pearson correlation for each Variable Methylated Region (VMR).
#'
#' This function will take a GRanges object converted into a data frame, where each row corresponds to a
#' Variable Methylated Region. Then, it computes the pairwise correlation of the probes of each VMR and reports
#' its median pairwise probe correlation.
#'
#' This function supports parallel computing for increased speed. To do so, you have to set the parallel backend
#' in your R session before running the function (e.g., doFuture::registerDoFuture()) and then the evaluation strategy (e.g., future::plan(multisession)). After that,
#' the function can be run as usual. It is recommended to also set options(future.globals.maxSize= +Inf).
#'
#' @inheritParams findVMRs
#' @return A data frame like VMR_df with an extra column per region containing the median pairwise correlation.
#'
#' @importFrom foreach %dopar%
#' @export
#'
medCorVMR = function(VMR_df, methylation_data){
  if(!is.list(VMRs_df$probes)){
    stop("Please make sure the 'probes' column in VMRs_df is a column of lists")
  }

  VMR_probes = VMR_df$probes #generate a list where each element will contain a vector with the probes present in one VMR
  #Compute correlations
  median_correlation = foreach::foreach(i = seq_along(VMR_probes), # For each VMR
                                        .combine = "c" #Combine outputs in a vector
                                        ) %dopar% {
    if (length(VMR_probes[[i]]) == 1){ #If the VMR has one probe
      NA
    }
    else{
      VMR_correlation = c()
      for (probe_x_i in 1:(length(VMR_probes[[i]])-1)){ #For each probe except the last one
        primary_probe = VMR_probes[[i]][probe_x_i]
        for (probe_y_i in (probe_x_i+1):length(VMR_probes[[i]])){ #compute the pairwise correlation with the downstream probes
          secondary_probe = VMR_probes[[i]][probe_y_i]
          VMR_correlation =  c(VMR_correlation,
                               stats::cor(unlist(methylation_data[primary_probe,]), #unlist added to make the subset df a vector
                                          unlist(methylation_data[secondary_probe,]),
                                          method= "pearson"))
        }
      }
      median_correlation = stats::median(VMR_correlation)
    }
  }

  VMR_df$median_correlation = median_correlation
  return(VMR_df)
}
