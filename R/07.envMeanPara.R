#' @title envMeanPara
#' @description
#' Calculate environmental mean parameters for a specific trait within a specified time range.
#' This time range is calculated by the Exhaustive_search.
#'
#' @param data Data frame containing trait mean values within different environments.
#' @param env_paras Data frame containing environmental parameter information.
#' @param maxR_dap1 The starting day for calculating environmental mean parameters.
#' Default is 18, which is used for the test data.
#' @param maxR_dap2 The ending day for calculating environmental mean parameters.
#' Default is 43, which is used for the test data.
#' @param Paras Vector containing the names of environmental parameters to be calculated.
#'
#' @return A data frame containing trait mean values and calculated environmental mean parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' envMeanPara <- envMeanPara(data = env_trait, env_paras = PTT_PTR,
#'                            maxR_dap1 = 18, maxR_dap2 = 43)
#' }
envMeanPara <- function( data, env_paras, maxR_dap1=NULL, maxR_dap2=NULL,Paras=NULL){

  # Set default values if not provided
  if (is.null(maxR_dap1)) {
    maxR_dap1 = 18
  }
  if (is.null(maxR_dap2)) {
    maxR_dap2 = 43
  }
  if(is.null(Paras)){
    Paras= c("DL", "GDD", "PTT", "PTR", "PTS")
  }

  # Get the number of environmental parameters
  nParas <- length(Paras)
  # Define the days range for calculation
  days <- c(maxR_dap1:maxR_dap2)

  # Initialize a matrix to store calculated environmental mean parameters
  env_facts_matrix <- matrix(nrow = nrow(data), ncol = nParas + 1)

  # Loop through each environment in the data
  for (e_i in 1:nrow(data)) {
    # Get the current environment code
    e <- data$env_code[e_i]

    # Subset environmental parameter data for the current environment
    env_para <- subset(env_paras, env_paras$env_code == e)

    # Calculate mean values for selected environmental parameters within the specified days
    env_mean <- colMeans(env_para[days, (1:nParas) + 4])  # Columns 5 and onwards contain parameter values

    # Assign calculated values to the matrix
    env_facts_matrix[e_i,] <- c(data$mean[e_i], round(env_mean, 4))
  }
  # Assign column names to the matrix
  colnames(env_facts_matrix) <- c('mean', Paras)

  # Merge calculated environmental mean parameters with the original data
  envMeanPara <- merge(data, env_facts_matrix)

  return(envMeanPara)
}
