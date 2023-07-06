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
#'
line_trait_mean <- function(data,LbyE,mean,trait,row){
  line_codes<-as.vector(unique(LbyE[,1]))
  lm_residuals <- data.frame(env_code = mean[,1])
  for (i in line_codes) {
    line_data_0 <- subset(data, data[["line_code"]] == i)# & !is.na(trait[["FTgdd"]]))
    if (nrow(line_data_0) > row) {
      line_data_0 <-  merge(env_trait, line_data_0)
      line_lm <- lm(line_data_0[[trait]] ~ line_data_0$mean)
      df1 <- data.frame(env_code = line_data_0$env_code, line_code = round((line_lm$residuals)^2, 3))
      colnames(df1)[2] <- i
      lm_residuals <- merge(lm_residuals, df1, all.x = T)
    }
  }
  data <- data.frame(env_code = lm_residuals[,1],
                     errors = rowMeans(lm_residuals[,-1], na.rm = T)) %>%
    left_join(env_trait)
  return(list(data,lm_residuals))
}

