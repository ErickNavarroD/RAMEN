#' Compute the median probe methylation correlation for each Variable Methylated Region
#'
#' This function will take a GRanges object converted into a data frame, where each row corresponds to a
#' Variable Methylated Region. Then, it computes the pairwise correlation of the probes of each VMR and report
#' its median pairwise probe correlation. OPTIONAL: This function supports parallel computing for increased speed. To do so, you have to set the parallel backend
#' in your R session before running the function (e.g., doFuture::registerDoFuture()) and then the evaluation strategy (e.g., future::plan(multisession)). After that,
#' the function can be run as usual.
#'
#' @param VMR_df GRanges object converted to a data frame. Must contain the following columns:
#' "seqnames", "start", "end"  (all of which are produced automatically when doing the object conversion)
#' and "probes" (containing a list where each element contains a vector with the probes
#' constituting the VMR).
#' @param data_methylation A data frame containing M or B values, with samples as columns and probes as rows
#'
#' @return A data frame like VMR_df with an extra column per region containing the median pairwise correlation
#'
#' @importFrom foreach %dopar%
#' @export
#'
medCorVMR = function(VMR_df, data_methylation){
  if(class(VMRs_df$probes) != "list"){
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
                               stats::cor(unlist(data_methylation[primary_probe,]), #unlist added to make the subset df a vector
                                          unlist(data_methylation[secondary_probe,]),
                                          method= "pearson"))
        }
      }
      median_correlation = median(VMR_correlation)
    }
  }

  VMR_df$median_correlation = median_correlation
  return(VMR_df)
}
