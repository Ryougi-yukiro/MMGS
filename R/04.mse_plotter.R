#' @title mse_plot
#' @description
#'
#' Generate a plot for Mean Squared Errors (MSE) data analyzed by line_trait_mean.
#'
#' @param data Dataframe containing MSE data analyzed by line_trait_mean.
#' @param point_size Size of points in the plot.
#' @param text_size Size of text labels in the plot.
#' @param point_shape Shape of the points in the plot.
#'
#' @return A plot illustrating Mean Squared Errors (MSE) data.
#' @export
#'
#' @importFrom ggplot2 aes geom_text geom_point element_blank
#' theme_bw theme unit
#'
#' @examples \dontrun{mse_plot(MSE)}
#' @usage mse_plot(data,point_size=NULL,text_size=NULL,point_shape=NULL)
mse_plot<-function(data,
                   point_size = NULL,
                   text_size = NULL,
                   point_shape = NULL){
  # Set default values if not provided
  if(is.null(point_size)){ point_size =3}
  if(is.null(text_size)){ text_size =3}
  if(is.null(point_shape)){ point_shape =20}
  ggplot(data)+
    geom_point(aes(data[,3], data[,2], color = data[,1]),
               size=point_size,shape=point_shape) +
    geom_text(aes(data[,3], data[,2], label = data[,1]),
              size = text_size,hjust=-1,vjust=1) +
    theme_bw(base_size = 8) +
    xlab("Mean")+
    ylab("Errors")+
    theme(panel.grid.major = element_blank(),
          legend.position = 'right',
          legend.key.size = unit(5, "pt"),
          panel.grid.minor = element_blank())
}
