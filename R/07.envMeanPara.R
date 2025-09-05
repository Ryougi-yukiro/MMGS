#' @title envMeanPara
#' @description
#' Calculate environmental mean parameters for a specific trait
#' within a specified time range.
#' This time range is calculated by the Exhaustive_search.
#'
#' @param data Data frame containing trait mean values
#' within different environments.
#' @param env_paras Data frame containing environmental parameter information.
#' @param max_d1 The starting day for calculating
#' environmental mean parameters.
#' Default is 18, which is used for the test data.
#' @param max_d2 The ending day for calculating
#' environmental mean parameters.
#' Default is 43, which is used for the test data.
#' @param Paras Vector containing the names of environmental
#' parameters to be calculated.
#'
#' @return A data frame containing trait mean values and calculated
#' environmental mean parameters.
#' @export
#'
#' @examples
#' \donttest{
#' env_trait<-env_trait_calculate(trait,"FTgdd","env_code")
#' envMeanPara<-EPM(data=env_trait,env_paras=PTT_PTR)
#' }
#'
EPM<-function(data,
           env_paras,
           max_d1=NULL,
           max_d2=NULL,
           Paras=NULL){
  # Set default values if not provided
  if (is.null(max_d1)) {
    max_d1 = 18
  }
  if (is.null(max_d2)) {
    max_d2 = 43
  }
  if(is.null(Paras)){
    Paras= c("DL", "GDD", "PTT", "PTR", "PTS")
  }

  # Get the number of environmental parameters
  nParas <- length(Paras)
  # Define the days range for calculation
  days <- c(max_d1:max_d2)

  # Initialize a matrix to store calculated environmental mean parameters
  env_facts_matrix <- matrix(nrow = nrow(data), ncol = nParas + 1)

  # Loop through each environment in the data
  for (e_i in 1:nrow(data)) {
    # Get the current environment code
    e <- data$env_code[e_i]

    # Subset environmental parameter data for the current environment
    env_para <- subset(env_paras, env_paras$env_code == e)

    env_mean <- colMeans(env_para[days, (1:nParas) + 4])

    # Assign calculated values to the matrix
    env_facts_matrix[e_i,] <- c(data$mean[e_i], round(env_mean, 4))
  }
  # Assign column names to the matrix
  colnames(env_facts_matrix) <- c('mean', Paras)

  # Merge calculated environmental mean parameters with the original data
  envMeanPara <- merge(data, env_facts_matrix)

  return(envMeanPara)
}
