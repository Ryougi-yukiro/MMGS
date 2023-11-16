#' Slope_Intercept
#'
#' @param maxR_dap1 as your like
#' @param maxR_dap2 as your like
#' @param Para_Name env para name
#' @param data trait data
#' @param line line code row
#' @param trait trait row
#' @param input env trait data
#' @param env_paras env para data
#' @param filter filter nums
#' @param rounds repeat nums
#'
#' @return Slope_Intercept dataframe
#' @export
#'
#' @examples \dontrun{Slope_Intercept<-Slope_Intercept(data=trait,input=env_trait, env_paras=PTT_PTR,
#'                                  Para_Name="PTT",line="line_code",trait="FTgdd",filter=5,
#'                                  maxR_dap1=18, maxR_dap2=43,rounds=4)}
Slope_Intercept <- function(maxR_dap1, maxR_dap2,Para_Name,data,line,trait,
                            input, env_paras,filter=NULL,rounds=NULL) {
  lines<-unique(as.vector(data[[line]]))
  PTT_PTR_ind <-  which(colnames(env_paras) == Para_Name);
  maxR_win <- c(maxR_dap1:maxR_dap2);
  prdM <- input;
  kPara <- c();
  for (e_i in 1:nrow(input)) {
    envParas <- subset(env_paras, env_paras$env_code == input$env_code[e_i]);
    envPara <- mean(envParas[maxR_win, PTT_PTR_ind]);
    kPara[e_i] <- envPara;
  }
  prdM$kPara <- kPara;
  lm_ab_matrix <- data.frame();
  for (l in 1: length(lines)) {
    l_trait <- subset(data, data[[line]] == lines[l]);
    if(nrow(l_trait) >= filter) {
      l_trait <- merge(l_trait, prdM);
      colnames(l_trait)[colnames(l_trait) == trait] <- "t_trait"
      lm_Mean <- lm(t_trait ~ mean, data = l_trait);
      lm_Para <- lm(t_trait ~ kPara, data = l_trait);
      a_Mean <- as.vector(round(predict(lm_Mean, data.frame(mean = mean(prdM$mean))), rounds)); ## adjusted by the population mean
      b_Mean <- as.vector(round(lm_Mean$coefficient[2], 4));
      b_Para <- as.vector(round(lm_Para$coefficient[2],4));
      a_Para <- as.vector(round(predict(lm_Para, data.frame(kPara = mean(prdM$kPara))), rounds)); ## adjusted by the population mean
      a_Para_ori <- as.vector(round(lm_Para$coefficient[1],4));
      merge <- c(lines[l], a_Mean, b_Mean, a_Para, b_Para);
      lm_ab_matrix<-rbind(lm_ab_matrix,merge)
    }

  }
  lm_ab_matrix <- lm_ab_matrix[!is.na(lm_ab_matrix[,2]),];
  colnames(lm_ab_matrix) <- c('line_code', 'Intcp_mean', 'Slope_mean',  'Intcp_para', 'Slope_para');
  return(lm_ab_matrix)
}
