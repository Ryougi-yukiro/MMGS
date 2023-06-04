#' ltm_plotter
#'
#' @param data line_trait_mean
#' @param mean env_trait
#' @param env_cols env colors gradient
#' @param shape point shape, like ggplot2
#' @param size point size, like ggplot2
#' @param linewidth line width , like ggplot2
#' @param area_color q25-q75 area fill color
#' @param area_alpha q25-q75 area fill color alpha
#'
#' @return A plot of population mean nad trait
#' @export
#'
#' @examples ltm_plotter(data=line_trait_mean,mean=env_trait)
ltm_plotter <- function(data,mean,env_cols =NULL,shape=NULL,
                        size=NULL,linewidth=NULL,area_color=NULL,area_alpha=NULL){
  if(is.null(env_cols)){
    env_cols <- rainbow_hcl(length(unique(as.vector(mean[,1]))), c = 80, l = 60,
                            start = 0, end = 300, fixup = TRUE, alpha = 0.75);
  }
  if(is.null(shape)){
    shape=18
  }
  if(is.null(size)){
    size=5
  }
  if(is.null(linewidth)){
    linewidth=1
  }
  if(is.null(area_color)){
    area_color= "#f0bfbb"
  }
  if(is.null(area_alpha)){
    area_alpha= 0.5
  }
  env_cols2<-factor(env_cols,levels = env_cols,labels = as.vector(mean[,1]))
  ggplot()+
    geom_line(data=data, aes(x = meanY , y = trait, group = line,
                             color = "gray80"),alpha=0.3,show.legend = FALSE)+
    geom_point(data=data, aes(x = meanY , y = trait, group = line,
                              color = "gray80"),alpha=0.1,show.legend = FALSE)+
    theme_bw()+labs(x = "mean", y = "trait")+
    geom_point(data = mean, aes(x = mean, y = mean,color = env_cols),
               size = size, shape = shape, show.legend = T)+
    geom_line(data = mean, aes(x = mean, y = mean, group = 1),
              color= "red",linewidth=linewidth,linetype='dashed', show.legend = FALSE)+
    geom_line(data = mean, aes(x = mean, y = q25, group = 1),
              color= "transparent",linewidth=1, show.legend = FALSE)+
    geom_line(data = mean, aes(x = mean, y = q75, group = 1),
              color= "transparent",linewidth=1, show.legend = FALSE)+
    geom_ribbon(data=mean,aes(x =mean, ymin = q25, ymax = q75,group = 1),
                fill = area_color,alpha=area_alpha)+
    theme(plot.title = element_text(size=12,hjust=0.5))+
    scale_color_manual(values = "gray80")
}

