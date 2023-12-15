#' @title env_trait_calculate
#' @description
#' Calculate the mean and optional quantiles for a trait across different environments.
#'
#' @param data Input data frame containing trait values and environmental information.
#' @param trait The name of the trait for which mean and optional quantiles will be calculated.
#' @param env The column name in the data frame representing different environments.
#' @param q25_75 Logical, indicating whether to calculate the 25th and 75th quantiles. Default is TRUE.
#'
#' @return A data frame containing the mean trait values for each environment. If q25_75 is TRUE, it also includes columns for the 25th (q25) and 75th (q75) quantiles, as well as the number of observations (n_obs).
#' @export
#'
#' @examples
#' # Calculate mean trait values and quantiles for "FTgdd" across different environments
#' env_trait <- env_trait_calculate(data = trait, trait = "FTgdd", env = "env_code", q25_75 = TRUE)
#'
#' # Access the result, e.g., mean trait values
#' env_trait$mean
#'
#' # Access the result, e.g., 25th quantile values
#' env_trait$q25
#'
#' # Access the result, e.g., 75th quantile values
#' env_trait$q75
#'
#' # Access the result, e.g., number of observations
#' env_trait$n_obs
env_trait_calculate<-function(data,trait,env,q25_75=TRUE){
  # Calculate mean trait values for each environment
  env_trait <- aggregate(x = data[[trait]],
                         by = list(env_code = data[[env]]),
                         mean)

  colnames(env_trait)[2] <- 'mean'

  if(q25_75 == TRUE){
    num_obs_env <- c();
    quantile_25 <- c();
    quantile_75 <- c();

    for (i in 1:nrow(env_trait)) {
     env_data <- subset(data, data[env] == env_trait[i,1])
     quantiles <- stats::quantile(env_data[[trait]], na.rm = T)
     quantile_25[i] <- quantiles[2]
     quantile_75[i] <- quantiles[4]
     num_obs_env[i] <- length(which(!is.na(env_data[,3])))
    }

    # Add quantiles and number of observations to the result
     env_trait$q25 <- quantile_25
     env_trait$q75 <- quantile_75
     env_trait$n_obs <- num_obs_env
  }
  return(env_trait)
}

