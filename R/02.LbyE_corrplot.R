#' @title LbyE_corrplot
#' @description
#' Plot the correlation between lines and environments.
#'
#' @param LbyE A data frame calculated by LbyE_calculate. Column names
#' represent Lines, and row names represent Environments.
#' @param cor_type Graphic format of the correlation map. The default format
#' is "heatmap". Other option is "pie".
#' @param color Set the color gradient for the correlation map. It is a
#' vector of three colors representing low, mid, and high values.
#'
#' @return A correlation plot between lines and environments in
#' the specified format.
#' @export
#' @importFrom
#' corrgram panel.shade panel.pie
#' @importFrom
#' ggplot2 aes geom_tile geom_text scale_fill_gradient2
#' theme_minimal xlab ylab ggplot
#' @examples
#' LbyE<-LbyE_calculate(trait,trait="FTgdd",env="env_code",line="line_code")
#' # Generate a heatmap correlation plot for Lines and Environments
#' LbyE_corrplot(LbyE,cor_type="heatmap",color=c("blue","white","red"))
#' # Generate a pie chart correlation plot for Lines and Environments
#' LbyE_corrplot(LbyE,cor_type="pie",color=c("blue","white","red"))

LbyE_corrplot <- function(LbyE, cor_type = NULL, color = NULL) {
  # Set default values if not provided
  if (is.null(cor_type)) {
    cor_type = "heatmap"
  }
  if (is.null(color)) {
    color = c("blue", "white", "red")
  }

  # Plot correlation based on the specified format
  if (cor_type == "pie") {
    corrgram::corrgram(as.matrix(LbyE[, -1]), order = TRUE,
                       pch = 19,
                       col.regions = grDevices::colorRampPalette(color))
  } else if (cor_type == "heatmap") {
    # Reshape data for ggplot
    corr_matrix <- stats::cor(as.matrix(LbyE[, -1]))

    # Reshape data for ggplot
    corr_data <- reshape2::melt(corr_matrix)

    # Create a ggplot for the heatmap
    p <- ggplot(corr_data, aes(x = corr_data[,1], y = corr_data[,2],
                               fill = corr_data[,3])) +
      geom_tile() +
      geom_text(aes(label = round(corr_data[,3], 2)), color = "white") +
      scale_fill_gradient2(low = color[1], mid = color[2], high = color[3], midpoint = 0) +
      theme_minimal() + xlab(NULL) + ylab(NULL)
    print(p)
  }
}
