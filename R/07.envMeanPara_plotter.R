#' @title envMeanPara_plotter
#' @description
#'
#' Generate a plot of environmental parameter means against the overall mean.
#'
#' @param data An environmental parameter mean data frame.
#' @param size Point size for scatter plot.
#' @param shape Point shape for scatter plot.
#' @param method Regression method for the smoothing line.
#' @param linewidth Line width for the smoothing line.
#' @param alpha Alpha value for points in the scatter plot.
#' @param linetype Line type for the smoothing line.
#' @param linecolor Line color for the smoothing line.
#' @param Paras Vector containing the names of environmental parameters to be considered.
#' @param env_code Point color for the scatter plot by env_code.Defualt is env_code.
#'
#' @return A plotter displaying the relationship between environmental parameter means and the overall mean.
#' @export
#' @importFrom ggplot2 geom_smooth annotate
#' @importFrom grDevices heat.colors
#'
#' @examples \dontrun{envMeanPara_plotter(envMeanPara,
#' Paras=c('DL', 'GDD', 'PTT', 'PTR', 'PTS'))}
envMeanPara_plotter <- function(data,size = NULL, shape = NULL,
                                method = NULL,Paras=NULL,env_code=NULL,
                                linewidth = NULL, alpha = NULL,
                                linetype = NULL, linecolor = NULL) {
  # Set default values if not provided
  if (is.null(size)) {
    size = 4
  }
  if (is.null(shape)) {
    shape = 20
  }
  if (is.null(method)) {
    method = "lm"
  }
  if (is.null(linewidth)) {
    linewidth = 1
  }
  if (is.null(alpha)) {
    alpha = 1
  }
  if (is.null(linetype)) {
    linetype = "dashed"
  }
  if (is.null(linecolor)) {
    linecolor = "black"
  }
  if (is.null(env_code)) {
    env_code = env_code
  }

  ps <- list()
  cols<-as.numeric(0)
  for (i in Paras) {
    df1<-data.frame()
    x <- (max(data[[i]]) + min(data[[i]])) / 2
    y <- max(data[, 1])
    r2 <- summary(lm(formula = data[[i]] ~ data[,1],
                     data = data))$r.squared
    df1 <- data.frame(x=data[,i],y=data[,1],env_code=data[,2])
    p <- ggplot() +
      geom_point(data=df1,
                 aes(x = x, y = y, color = env_code),
                 size = size, shape = shape, alpha = alpha) +
      geom_smooth(data=df1,
                  aes(x = x, y = y),
                  se = FALSE, linewidth = linewidth,
                  method = method, linetype = linetype, color = linecolor, formula = y ~ x) +
      xlab(colnames(data)[cols])+ylab("Mean")+xlab(i)+
      theme_bw() +
      annotate('text', x = x, y = y,
               label = (paste("R2 =", round(r2, 2))),
               size = 3, vjust = 0.5, hjust = 0.5)

    ps[[i]] <- p
  }
  d <- character(length(ps))

  for (i in 1:length(ps)) {
    d[i] <- paste0("ps[[", i, "]]")
  }

  a <- paste(d, collapse = ", ")
  code <- paste("lemon::grid_arrange_shared_legend(", paste(a, collapse = ", "), ", nrow = 2, ncol = 3, position = 'right')")
  eval(parse(text = code))
}
