#' @title etl_plotter
#' @description
#' Generate a plot of environmental trait lines with mean, q25, and q75 areas.
#'
#' @param data An etl input data frame.
#' @param trait An env_trait data frame containing mean, q25, and q75 values.
#' @param env_cols Environmental colors gradient; if NULL, it generates colors using rainbow_hcl.
#' @param shape Point shape, similar to ggplot2.
#' @param size Point size, similar to ggplot2.
#' @param linewidth Line width, similar to ggplot2.
#' @param area_color Fill color for the area between q25 and q75.
#' @param area_alpha Fill color alpha for the area between q25 and q75.
#'
#' @return A plot of etl, filling the area between env-trait q25 to q75, with a red line representing the env-trait mean.
#' @export
#'
#' @importFrom colorspace rainbow_hcl
#' @importFrom ggplot2 element_text geom_line geom_point geom_ribbon scale_color_manual
#'
#' @examples
#' #Get input:
#' env_trait<-env_trait_calculate(data=trait,trait="FTgdd",env="env_code")
#' LbyE<-LbyE_calculate(data=trait,trait="FTgdd",env="env_code",line="line_code")
#' etl<-LbyE_Reshape(data=env_trait,LbyE=LbyE,env="env_code")
#' #Plot
#' etl_plotter(data=etl,trait=env_trait)
etl_plotter <- function(data, trait, env_cols = NULL, shape = NULL,
                        size = NULL, linewidth = NULL, area_color = NULL, area_alpha = NULL) {
  # Set default values if not provided
  if (is.null(env_cols)) {
    env_cols <- colorspace::rainbow_hcl(length(unique(as.vector(data[, 1]))),
                                        c = 80, l = 60,
                                        start = 0, end = 300,
                                        fixup = TRUE, alpha = 0.75)
  }
  if (is.null(shape)) {
    shape <- 18
  }
  if (is.null(size)) {
    size <- 5
  }
  if (is.null(linewidth)) {
    linewidth <- 1
  }
  if (is.null(area_color)) {
    area_color <- "#f0bfbb"
  }
  if (is.null(area_alpha)) {
    area_alpha <- 0.5
  }

  p<-ggplot() +
    geom_line(data = data, aes(x = data[,1], y = data[,2] , group = data[,3],
                               color = "gray80", alpha = 0.1), show.legend = FALSE) +
    geom_point(data = data, aes(x = data[,1], y = data[,2] , group = data[,3],
                                color = "gray80", alpha = 0.5), alpha = 0.1, show.legend = FALSE) +
    theme_bw()
  p <-p + labs(x = "env", y = "Trait") +
    geom_point(data = trait, aes(x = trait[,1], y = trait[,2]),
               color = env_cols, size = size, shape = shape, show.legend = FALSE) +
    geom_point(data = trait, aes(x = trait[,1], y = trait[,3]),
               color = env_cols, size = size, shape = shape, show.legend = FALSE) +
    geom_point(data = trait, aes(x = trait[,1], y = trait[,4]),
               color = env_cols, size = size, shape = shape, show.legend = FALSE) +
    geom_line(data = trait, aes(x = trait[,1], y = trait[,2], group = 1),
              color = "red", linewidth = linewidth, show.legend = FALSE) +
    geom_line(data = trait, aes(x = trait[,1], y = trait[,3], group = 1),
              color = "black", linewidth = linewidth, show.legend = FALSE) +
    geom_line(data = trait, aes(x = trait[,1], y = trait[,4], group = 1),
              color = "black", linewidth = linewidth, show.legend = FALSE) +
    geom_ribbon(data = trait, aes(x = trait[,1], ymin = trait[,3], ymax = trait[,4], group = 1),
                fill = area_color, alpha = area_alpha) +
    theme(plot.title = element_text(size = 12, hjust = 0.5)) +
    scale_color_manual(values = "gray80")
  return(p)
}
