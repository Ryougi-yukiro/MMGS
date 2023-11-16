#' envMeanPara_plotter
#'
#' @param envMeanPara  envMeanPara
#' @param size point size
#' @param shape point shape
#' @param method regression method
#' @param linewidth regression line width
#' @param alpha point alpha
#' @param linetype regression line type
#' @param linecolor regression line color
#'
#' @return A envparamean plotter
#' @export
#'
#' @examples \dontrun{envMeanPara_plotter(envMeanPara)}
envMeanPara_plotter<-function(envMeanPara,size=NULL,shape=NULL,method=NULL,
                              linewidth=NULL,alpha=NULL,linetype=NULL,linecolor=NULL){
  if(is.null(size)){
    size= 4
  }
  if(is.null(shape)){
    shape= 20
  }
  if(is.null(method)){
    method= "lm"
  }
  if(is.null(linewidth)){
    linewidth= 1
  }
  if(is.null(alpha)){
    alpha= 1
  }
  if(is.null(linetype)){
    linetype= "dashed"
  }
  if(is.null(linecolor)){
    linecolor= "black"
  }
  ps <- list()
  for (i in Paras){
    x<-(max(envMeanPara[[i]])+min(envMeanPara[[i]]))/2;
    y<-max(envMeanPara[,1]);
    r2<-summary(lm(formula = envMeanPara[[i]] ~ mean, data = envMeanPara))$r.squared
    p<-ggplot()+geom_point(data=envMeanPara,aes(x=.data[[i]],y=mean,color=env_code),
                            size=size,shape=shape,alpha=alpha)+
      geom_smooth(data=envMeanPara,aes(x=.data[[i]],y=mean),se=F,linewidth=linewidth,method=method,
                  linetype=linetype,color=linecolor,formula = y ~ x)+
      theme_bw()+
      annotate('text',x=x,y=y,
              label = (paste("R2 =", round(r2, 2))),
              size=3,vjust = 0.5,hjust = 0.5)

    ps[[i]] <- p
  }
  d <- character(length(ps))


  for (i in 1:length(ps)) {d[i] <- paste0("ps[[", i,"]]")}

  a<-paste(d, collapse = ", ")
  code <- paste("lemon::grid_arrange_shared_legend(", paste(a, collapse = ", "), ", nrow = 2, ncol = 3, position = 'right')")
  eval(parse(text = code))

}
