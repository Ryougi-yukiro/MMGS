#' etl_plotter
#'
#' @param data etl input
#' @param mean env_trait
#' @param env_cols env colors gradient
#' @param shape point shape, like ggplot2
#' @param size point size, like ggplot2
#' @param linewidth line width , like ggplot2
#' @param area_color q25-q75 area fill color
#' @param area_alpha q25-q75 area fill color alpha
#'
#' @return A plot of etl ,fill area start env-trait q25 to q75, red line is env-trait mean
#' @export
#'
#' @examples etl_plotter(data=etl,mean=env_trait)
etl_plotter <- function(data,mean,env_cols =NULL,shape=NULL,
                        size=NULL,linewidth=NULL,area_color=NULL,area_alpha=NULL){
  if(is.null(env_cols)){
    env_cols <- colorspace::rainbow_hcl(length(unique(as.vector(data[,1]))), c = 80, l = 60,
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
  ggplot()+
    geom_line(data=data, aes(x = env_code , y = trait, group = line,
                             color = "gray80",alpha=0.1),show.legend = FALSE)+
    geom_point(data=data, aes(x = env_code , y = trait, group = line,
                              color = "gray80",alpha=0.5),alpha=0.1,show.legend = FALSE)+
    theme_bw()+labs(x = "env", y = "trait")+
    geom_point(data = mean, aes(x = env_code, y = mean),
               color = env_cols, size = size, shape = shape, show.legend = FALSE)+
    geom_point(data = mean, aes(x = env_code, y = q25),
               color = env_cols, size = size, shape = shape, show.legend = FALSE)+
    geom_point(data = mean, aes(x = env_code, y = q75),
               color = env_cols, size = size, shape = shape, show.legend = FALSE)+
    geom_line(data = mean, aes(x = env_code, y = mean, group = 1),
              color= "red",linewidth=linewidth, show.legend = FALSE)+
    geom_line(data = mean, aes(x = env_code, y = q25, group = 1),
              color= "black",linewidth=linewidth, show.legend = FALSE)+
    geom_line(data = mean, aes(x = env_code, y = q75, group = 1),
              color= "black",linewidth=linewidth, show.legend = FALSE)+
    geom_ribbon(data=mean,aes(x =env_code, ymin = q25, ymax = q75,group = 1),
                fill = area_color,alpha=area_alpha)+
    theme(plot.title = element_text(size=12,hjust=0.5))+
    scale_color_manual(values = "gray80")
}
