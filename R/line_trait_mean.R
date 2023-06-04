#' line_trait_mean
#'
#' @param LbyE LbyE input
#' @param mean env_trait input
#' @param row mean rows in env_trait
#'
#' @return line_trait_mean
#' @export
#'
#' @examples line_trait_mean(LbyE=LbyE,mean=env_trait,row=2)
line_trait_mean <- function(LbyE,mean,row){

data<-data.frame()
for (i in 1:nrow(LbyE)) {
  df <- data.frame(meanY = mean[,row], trait = as.numeric(LbyE[i, -1]));
  df <- df[!is.na(df$trait),];
  df$line<-LbyE[i, 1]
  data<-rbind(data,df)
}
return(data)
}
