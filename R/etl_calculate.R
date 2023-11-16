#' etl_calculate
#' @description
#' calculate your trait values within one env and one line
#'
#' @param data env info data
#' @param trait trait you select
#' @param env your envs col
#' @param line your lines col
#' @param bycol order your env data by your select.
#'
#' @return your trait values within one env and one line
#' @export
#'
#' @examples \dontrun{etl<-etl_calculate(data=env_info,trait="FTgdd",env="env_code",bycol="lat")}

etl_calculate <- function(data,trait,env,line,bycol){
  env_trait <- env_trait[env_trait[[env]]%in% colnames(LbyE),]

  env_geo <- merge(env_trait, env_info);

  env_geo <- env_geo[order(env_geo[[bycol]]),];

  env_geo_order <- match(env_trait[[env]], env_trait[[env]]);

  data<-data.frame()
  for (i in 1:nrow(LbyE)) {
    df <- data.frame(env_code = env_trait[[env]],
                      env_order = env_geo_order,trait = as.numeric(LbyE[i, -1]));
    df<- df[!is.na(df$trait),];
    df <- df[order(df$env_order),];
    df$line <- LbyE[i,1]
    data<-rbind(data,df)
  }
  return(data)
}
