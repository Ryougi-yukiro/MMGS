#' @title LbyE_calculate
#' @description
#' calculate Lines in all envs
#'
#' @param env All envs need calculate
#' @param line All Lines need calculate
#' @param data Input data
#' @param trait Trait you select
#'
#' @return Lines values within Envs.
#' @export
#'
#' @examples \dontrun{LbyE<-LbyE_calculate(data=data,trait="FTgdd",env="env_code",line="line_code")}
LbyE_calculate<-function(data,trait,env,line){
  #输出每个line在每个环境中平均表型
  envs <- unique(as.vector(data[[env]]));
  lines <- unique(as.vector(data[[line]]));
  LbyE <- data.frame(line_code = lines);
  t<-match(trait, colnames(data))
  t2<-match(line, colnames(data))
  for (e_l in 1:length(envs)) {
    e <- envs[e_l];
    e_l_trait <- subset(data, data[[env]] == e);
    nonNAs <- length(which(!is.na(e_l_trait[,t])))
    colnames(e_l_trait)[t] <- e;
    LbyE <- merge(LbyE, e_l_trait[,c(t2,t)], all.x = T)
    #LbyE <- merge(LbyE, e_trait[,c(1,match(env, colnames(e_trait)))], all.x = T)
  }
  return(LbyE)
}
