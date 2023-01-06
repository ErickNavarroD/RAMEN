#' Compute the median correlation of a Genomic Ranges data frame
#'
#' Given a set of candidate VMRs, this function will compute the pairwise correlation of the probes
#'  in each VMR and report its median.
#'
#' @param VMR_df A data frame converted from a Genomic Ranges object. Must contain the following columns: "probes", containing a list where each
#' element contains a vector with the probes constituting a  VMR
#' @param data_methylation A data frame containing M or B values, with samples as columns and probes as rows
#'
#' @return A data frame like VMR_df, but with an extra column per region containing the median pairwise correlation
#' @export
#'
medCorVMR = function(VMR_df, data_methylation){
  median_correlation = c() #Start the vector that will be outputted
  probes = VMR_df$probes #generate a list where each element will contain a vector with the probes present in one VMR
  #Compute correlations
  for (VMR in probes){
    if (length(VMR) == 1){
      median_correlation = c(median_correlation, NA)
    }
    else{
      VMR_correlation = c()
      for (probe_x_i in 1:(length(VMR)-1)){ #For each probe except the last one
        primary_probe = VMR[probe_x_i]
        for (probe_y_i in (probe_x_i+1):length(VMR)){ #compute the pairwise correlation with the downstream probes
          secondary_probe = VMR[probe_y_i]
          VMR_correlation =  c(VMR_correlation,
                               stats::cor(unlist(data_methylation[primary_probe,]), #unlist added to make the subset df a vector
                                          unlist(data_methylation[secondary_probe,]),
                                          method= "pearson"))
        }
      }
      median_correlation = c(median_correlation, stats::median(VMR_correlation))
    }
  }
  VMR_df$median_correlation = median_correlation
  return(VMR_df)
}
