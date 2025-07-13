#' Mean_trait_plot
#'
#' Generate a plot displaying the relationship between the mean trait values
#' and the environmental parameter mean.
#'
#' @param Reg Regression result data frame.
#' @param MSE Mean Squared Error (MSE) result data frame.
#' @param point_size Size of points in the plot.
#' @param point_shape Shape of points in the plot.
#' @param linewidth Line width for the smoothing line.
#' @param linetype Line type for the smoothing line.
#' @param color color for the points and the dash area.
#'
#' @return A plot showing the relationship between mean trait values and
#' the environmental parameter mean, with MSE information.
#'
#' @export
#'
#' @importFrom dplyr %>% left_join
#' @importFrom ggplot2 geom_line geom_ribbon geom_point geom_smooth theme_bw
#' @examples \dontrun{Mean_trait_plot(Reg, MSE)}
#'
Mean_trait_plot <- function(Reg, MSE, point_size = NULL,color=NULL,
                            point_shape = NULL, linewidth = NULL,
                            linetype = NULL) {
  # Set default values if not provided
  if (is.null(point_size)) {
    point_size = 0.5
  }
  if (is.null(point_shape)) {
    point_shape = 1
  }
  if (is.null(linewidth)) {
    linewidth = 0.5
  }
  if (is.null(linetype)) {
    linetype = 'dashed'
  }
  if (is.null(color)) {
    color = c("grey","#EE82EE37")
  }
  data <- Reg %>% left_join(MSE)
  ggplot(data) +
    geom_line(aes(data[,1], data[,2], group = data[,4]),
              color = "#80808023", show.legend = FALSE) +
    geom_point(aes(data[,1], data[,2], group = data[,4]),
               color = "#80808023", show.legend = FALSE) +
    geom_point(aes(data[,1], data[,2], group = data[,4],
                   color = data[,3]), size = point_size, shape = point_shape) +
    geom_smooth(aes(data[,1], data[,2]), color = color[1],
                linewidth = linewidth,
                method = 'lm', se = FALSE, linetype = linetype) +
    geom_ribbon(aes(data[,1], ymin = data[,7], ymax = data[,8]),
                fill = color[2]) +
    geom_point(aes(data[,1], data[,1])) + theme_bw()
}
