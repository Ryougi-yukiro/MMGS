#' LOOCV
#'
#' @param data trait data
#' @param input env trait data
#' @param env_paras env para info
#' @param Para_Name Para name
#' @param line line name
#' @param trait trait name
#' @param maxR_dap1 as you like depend on your data
#' @param maxR_dap2 as you like depend on your data
#' @param p defualt=1
#' @param filter filter num
#'
#' @return LOOCV data frame
#' @export
#'
#' @examples \dontrun{LOOCV(maxR_dap1=18,maxR_dap2=43,data=trait,input=env_trait,env_paras=PTT_PTR,
#'                Para_Name="PTT",trait="FTgdd",line="line_code",p=1,filter=4)}
LOOCV<-function(data,input,env_paras,Para_Name,line,
                trait,maxR_dap1=NULL,maxR_dap2=NULL,p=NULL,filter=NULL){
  maxR_win <- c(maxR_dap1:maxR_dap2);
  PTT_PTR_ind <-  which(colnames(env_paras) == Para_Name);
  prdM <- env_trait;
  maxR_envPara <- matrix(ncol = 2, nrow = nrow(env_trait));
  kPara <- c();
  for (e_i in 1:nrow(env_trait)) {
    envParas <- subset(env_paras, env_paras[,1] == env_trait[,1][e_i]);
    envPara <- mean(envParas[maxR_win, PTT_PTR_ind]);
    kPara[e_i] <- envPara;
  }
  prdM$kPara <- kPara;
  obs_prd_m <- matrix(0, ncol = 7,nrow = nrow(data));
  n <- 0;
  lines<-unique(as.vector(data[[line]]))
  for (l in lines) {
    l_trait <- subset(data, data[[line]] == l);
    ril_data <- merge(prdM, l_trait,  all.x = T);
    colnames(ril_data)[colnames(ril_data) == trait] <- "t_trait"
    if (length(which(!is.na(ril_data[["t_trait"]]))) > 4) {
      for (e in 1:nrow(ril_data)) {
        obs_trait <- ril_data[["t_trait"]][e];
        if (!is.na(obs_trait)) {
          trn <- ril_data[-e,];
          l_mean <- mean(trn[["t_trait"]], na.rm = T);
          prd_trait_mean  <- round(predict( lm(t_trait ~ mean, data = trn), ril_data[e,]), digits = 3);
          prd_trait_kpara <- round(predict( lm(t_trait ~ kPara, data = trn), ril_data[e,]), digits = 3);
          n <- n + 1;
          obs_prd_m[n,] <- c(ril_data[,1][e], p, l,
                             prd_trait_mean, prd_trait_kpara, obs_trait, l_mean);
        }
      }
    }
  }
  obs_prd_m <- obs_prd_m[1:n,]
  colnames(obs_prd_m) <- c('env_code', 'pop_code', 'ril_code',
                           'Prd_trait_mean', 'Prd_trait_kPara', 'Obs_trait', 'Line_mean');
  obs_prd_m<-as.data.frame(obs_prd_m)
  return(obs_prd_m);
}
