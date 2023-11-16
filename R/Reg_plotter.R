#' Title
#'
#' @param Reg input data
#' @param size piont size
#' @param method smooth method
#' @param color line and point color
#' @param alpha line and point alpha
#'
#' @return A regression plot between mean and trait
#' @export
#' @importFrom ggplot2 aes geom_smooth geom_point scale_fill_gradient2 labs theme_bw
#'
#' @examples \dontrun{Reg_plotter(Reg=Reg)}
Reg_plotter<-function(Reg=Reg,size=NULL,method=NULL,color=NULL,alpha=NULL){
  if(is.null(size)){
    size= 1
  }
  if(is.null(method)){
    method= "lm"
  }
  if(is.null(color)){
    color= c("gray90","gray50")
  }
  if(is.null(alpha)){
    alpha= c(0.1,0.5)
  }
  ggplot(Reg, aes(x = meanY, y = Trait, group = line)) +
    geom_smooth(method = method, se = FALSE,color=color[1],alpha=alpha[1]) +
    geom_point(color=color[2],size=size,alpha=alpha[2])+labs(x = "mean", y = "Trait") +theme_bw()


}
