#' Regression
#'
#' @param LbyE LbyE
#' @param env_trait env trait info
#' @param filter_num filter num
#'
#' @return Regression calculate result
#' @export
#'
#' @examples Reg(LbyE=LbyE,env_trait=env_trait)
Reg<-function(LbyE,env_trait,filter_num=NULL){
  if(is.null(filter_num)){
    filter_num= 4
  }
  data<-data.frame()
  for (i in 1:nrow(LbyE)) {
    df <- data.frame(meanY = env_trait$mean, Trait = as.numeric(LbyE[i, -1]));
    df <- df[!is.na(df$Trait),];
    if(nrow(df) >= filter_num) {
      df$line<-LbyE[i,1]
      data<-rbind(data,df)
    }
  }
  return(data)
}
