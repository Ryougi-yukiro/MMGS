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
<<<<<<< HEAD
#' @examples \dontrun{ltm_plotter(data=line_trait_mean,mean=env_trait)}
=======
#' @examples ltm_plotter(data=line_trait_mean,mean=env_trait)
>>>>>>> 23b0efae316a77181b94cb5eff39d6257393dcf5
ltm_plotter <- function(data,mean,env_cols =NULL,shape=NULL,
                        size=NULL,linewidth=NULL,area_color=NULL,area_alpha=NULL){
  if(is.null(env_cols)){
    env_cols <- colorspace::rainbow_hcl(length(unique(as.vector(mean[,1]))), c = 80, l = 60,
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
  ggplot(Reg %>% left_join(ltm)) +
    # geom_line(aes(env_order, Yobs, group = g), color = "#80808023", show.legend = F) +
    geom_line(aes(meanY, Yobs, g = g), color = "#80808023", show.legend = F) +
    geom_point(aes(meanY, Yobs, group = g), color = "#80808023", show.legend = F) +
    geom_point(aes(meanY, Yobs, group = g, color = env_code),size = 0.5, shape = 1) +
    geom_smooth(aes(meanY, Yobs), color = 'grey', size = 0.5, method = 'lm', se = F, linetype = 'dashed') +
    geom_ribbon(aes(meanY, ymin = q1, ymax = q3), fill = "#EE82EE37") +
    geom_point(aes(meanY, meanY)) +
    theme_main +
    xlab('Population mean') +
    ylab('FTgdd')
}

