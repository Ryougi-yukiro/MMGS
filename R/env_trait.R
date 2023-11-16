#' @title env_trait_calculate
#' @description
#' Calculate traits mean in all envs
#'
#' @param data Input data
#' @param trait The trait need calculate
#' @param env All envs need calculate
#' @param q25_75 Calculate q25 and q75 values. The defualt is TRUE.
#'
#' @return Trait mean values within different envs. if q25_75=T, also contains q25 and q75 values
#' @export
#'
#' @examples \dontrun{env_trait<-env_trait(data=trait,trait="FTgdd",env="env_code")}
env_trait_calculate<-function(data,trait,env,q25_75=TRUE){
  env_trait <- aggregate(x = data[[trait]],by = list(env_code = data[[env]]),mean);
  colnames(env_trait)[2] <- 'mean';

  env_trait <- env_trait[order(env_trait$mean),];

  if(q25_75 == TRUE){
    num_obs_env <- c();
    quantile_25 <- c();
    quantile_75 <- c();

    for (i in 1:nrow(env_trait)) {
     env_data <- subset(data, data[env] == env_trait[i,1]);
     quantiles <- stats::quantile(env_data[[trait]], na.rm = T);
     quantile_25[i] <- quantiles[2];
     quantile_75[i] <- quantiles[4];
     num_obs_env[i] <- length(which(!is.na(env_data[,3])))
    }
     env_trait$q25 <- quantile_25;
     env_trait$q75 <- quantile_75;
     env_trait$n_obs <- num_obs_env;
  }
  return(env_trait)
}

