#' Exhaustive_search
#'
#' @param data env_trait data from env_trait_calculate
#' @param env_paras env_paras data
#' @param searching_daps depends on your file,set your para, default is 122 depends on example data
#' @param p defualt is 1
#' @param dap_x like searching_daps
#' @param dap_y like searching_daps
#' @param LOO defualt is 0
#' @param Paras your env
#'
#' @return pop cor matrix
#' @export
#'
#' @examples \dontrun{pop_cor<-Exhaustive_search(data=env_trait, env_paras=PTT_PTR, searching_daps=122,
#'                                      p=1, dap_x=122,dap_y=122,LOO=0,Paras=Paras)}
Exhaustive_search<-function(data, env_paras, searching_daps=NULL,
                            p=NULL, dap_x=NULL,dap_y=NULL,LOO=NULL,Paras){
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
  nParas <- length(Paras);
  dap_win <- searching_daps * searching_daps  / 2;
  pop_cors_matrix <- matrix(ncol = 4 + (2 * nParas), nrow = dap_win * 1);
  colnames(pop_cors_matrix) <- c("pop_code", 'Day_x', 'Day_y', 'window',
                                 paste('R_', Paras, sep = ''), paste('nR_', Paras, sep = ''));
  n <- 0;
  for (d1 in 1:(dap_y - 6)) {
    for (d2 in (d1 + 6):dap_y) {
      days <- c(d1:d2);
      env_facts_matrix <- matrix(nrow = nrow(data), ncol = nParas);
      for (e_i in 1:nrow(data)) {
        e <- data$env_code[e_i];
        env_para <- subset(env_paras, env_paras$env_code == e);
        env_mean <- colMeans(env_para[days, (1:nParas) + 4]);
        env_facts_matrix[e_i,] <- env_mean;

      }
      n <- n + 1;
      Ymean_envPara <- cbind(env_facts_matrix, data$mean);
      #print(Ymean_envPara)
      rs <- c();
      if (LOO == 0) {
        for (k in 1:nParas) {
          rs[k] <- round(cor(Ymean_envPara[,nParas + 1], Ymean_envPara[,k]), digits = 4)

        }
      } else {
        loo_rs_matrix <- matrix(nrow = nrow(Ymean_envPara)+ 0, ncol = nParas);
        for (k in 1:nParas) { ## 8 environment parameters
          for (e_x in c(1:nrow(Ymean_envPara))) {
            t_matrix <- Ymean_envPara[-e_x,];
            loo_rs_matrix[e_x, k] <- round(cor(t_matrix[,nParas + 1], t_matrix[,k]), digits = 4)
          }
        }
        rs <- apply(loo_rs_matrix, 2, median);
      }
      pop_cors_matrix[n, ] <- c(p, d1, d2, d2 - d1, rs, 0 - rs);
    }
  }
  pop_cors_matrix <- pop_cors_matrix[1:n,]
  return(pop_cors_matrix)
}
