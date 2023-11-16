#' mse_plotter
#'
#' @param data mse data analysised by line_trait_mean
#' @param point_size point size
#' @param text_size text size
#' @param point_shape point shape
#'
#' @return A ltm ploter
#' @export
#' @importFrom ggplot2 aes geom_text geom_point element_blank theme_bw theme
#' @examples
#' \dontrun{mse_plotter(ltm)}
mse_plotter<-function(data,point_size=NULL,text_size=NULL,point_shape=NULL){
  if(is.null(point_size)){ point_size =3}
  if(is.null(text_size)){ text_size =3}
  if(is.null(point_shape)){ point_shape =20}
  ggplot(data)+
  geom_point(aes(mean, errors, color = env_code),size=point_size,shape=point_shape) +
  geom_text(aes(mean, errors, label = env_code), size = text_size) +
  theme_bw(base_size = 8) +
  theme(panel.grid.major = element_blank(),
        legend.position = 'right',
        legend.key.size = unit(5, "pt"),
        panel.grid.minor = element_blank())
}
