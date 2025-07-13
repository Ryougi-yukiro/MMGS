#' @title Exhaustive_search
#' @description
#' Perform exhaustive search for population
#' correlation matrix based on environmental parameters.
#'
#' @param data Data frame containing environmental trait data.
#' @param env_paras Data frame containing environmental parameter information.
#' @param searching_daps The total number of days to search for,
#' typically set based on your data. Default is 122.
#' @param p The parameter for controlling the shape of the search window.
#' Default is 1.
#' @param dap_x The number of days for the x-axis of the search window.
#'  Default is the same as searching_daps.
#' @param dap_y The number of days for the y-axis of the search window.
#' Default is the same as searching_daps.
#' @param LOO Leave-One-Out cross-validation flag. If LOO is 1,
#  it performs Leave-One-Out cross-validation. Default is 0.
#' @param Paras Vector containing the names of environmental parameters
#' to be considered.
#'
#' @return A population correlation matrix based on the exhaustive search.
#' @export
#'
#' @examples
#' #Input from function : env_trait_calculate
#' env_trait<-env_trait_calculate(data=trait,trait="FTgdd",env="env_code")
#' #Run
#' \dontrun{Exhaustive_search(denv_trait,PTT_PTR)}
#'
#' @usage Exhaustive_search(data, env_paras, searching_daps, ...)
Exhaustive_search<-function(data, env_paras, searching_daps,
                            p=NULL, dap_x=NULL,dap_y=NULL,LOO=NULL,
                            Paras=NULL){
  if(is.null(p)){
    p= 1
  }
  if(is.null(dap_x)){
    dap_x= searching_daps
  }
  if(is.null(dap_y)){
    dap_y= searching_daps
  }
  if(is.null(LOO)){
    LOO= 0
  }
  if(is.null(Paras)){
    Paras= c("DL", "GDD", "PTT", "PTR", "PTS")
  }
  nParas <- length(Paras);
  dap_win <- searching_daps * searching_daps / 2;
  pop_cors_matrix <- matrix(ncol = 4 + (2 * nParas), nrow = dap_win * 1);
  colnames(pop_cors_matrix) <- c("pop_code", 'Day_x', 'Day_y', 'window',
                                 paste('R_', Paras, sep = ''),
                                 paste('nR_', Paras, sep = ''));
  n <- 0;

  # Loop through each window of days
  for (d1 in 1:(dap_y - 6)) {
    for (d2 in (d1 + 6):dap_y) {
      days <- c(d1:d2);
      env_facts_matrix <- matrix(nrow = nrow(data), ncol = nParas);

      # Calculate environmental means for each environment
      for (e_i in 1:nrow(data)) {
        e <- data$env_code[e_i];
        env_para <- subset(env_paras, env_paras$env_code == e);
        env_mean <- colMeans(env_para[days, (1:nParas) + 4]);
        env_facts_matrix[e_i,] <- env_mean;
      }

      n <- n + 1;
      Ymean_envPara <- cbind(env_facts_matrix, data$mean);
      rs <- c();

      # Calculate correlations
      if (LOO == 0) {
        for (k in 1:nParas) {
          rs[k] <- round(cor(Ymean_envPara[, nParas + 1],
                             Ymean_envPara[, k]), digits = 4)
        }
      } else {
        loo_rs_matrix <- matrix(nrow = nrow(Ymean_envPara) + 0, ncol = nParas);
        for (k in 1:nParas) {
          ## 8 environment parameters
          for (e_x in c(1:nrow(Ymean_envPara))) {
            t_matrix <- Ymean_envPara[-e_x,];
            loo_rs_matrix[e_x, k] <- round(cor(t_matrix[, nParas + 1],
                                               t_matrix[, k]), digits = 4)
          }
        }
        rs <- apply(loo_rs_matrix, 2, median);
      }

      # Populate the population correlation matrix
      pop_cors_matrix[n, ] <- c(p, d1, d2, d2 - d1, rs, 0 - rs);
    }
  }

  # Trim the matrix to the actual number of windows
  pop_cors_matrix <- pop_cors_matrix[1:n,]
  return(pop_cors_matrix)
}
