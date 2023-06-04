#' envMeanPara
#'
#' @param data env_trait data
#' @param env_paras env para info data
#' @param maxR_dap1 as your like
#' @param maxR_dap2 as your like
#' @param Paras All env para names
#'
#' @return A dat frame of envMeanPara
#' @export
#'
#' @examples envMeanPara<-envMeanPara(data=env_trait, env_paras=PTT_PTR, maxR_dap1=18,maxR_dap2=43, Paras=Paras)
envMeanPara <- function( data, env_paras, maxR_dap1=NULL, maxR_dap2=NULL,Paras){
  if(is.null(maxR_dap2)){
    maxR_dap1= 18
  }
  if(is.null(maxR_dap2)){
    maxR_dap1= 43
  }
  nParas <- length(Paras);
  days <- c(maxR_dap1:maxR_dap2);
  env_facts_matrix <- matrix(nrow = nrow(data), ncol = nParas + 1);
  for (e_i in 1:nrow(data)) {
    e <- data$env_code[e_i];
    env_para <- subset(env_paras, env_paras$env_code == e);
    env_mean <- colMeans(env_para[days, (1:nParas) + 4]); ### DL, GDD, PTT, PTR, PTS

    env_facts_matrix[e_i,] <- c(data$mean[e_i], round(env_mean, 4) );
  }
  colnames(env_facts_matrix) <- c( 'mean', Paras);
  envMeanPara <- merge(data, env_facts_matrix);
  return(envMeanPara)
}
