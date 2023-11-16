#' @title LbyE_corrplot
#' @description
#' Plot Lines and envs correlation
#'
#' @param LbyE A data frame calculate by LbyE_calculate, the colnames are Lines, the rawnames are Envs.
#' @param cor_type Graphic format of the corr map.The defualt format is heatmap.
#' @param color Set the color gradient
#'
#' @return A correalation plot between lines and envs
#' @export
#' @importFrom corrgram panel.shade panel.pie
#' @importFrom ggplot2 aes geom_tile geom_text scale_fill_gradient2 theme_minimal xlab ylab ggplot
#' @examples \dontrun{LbyE_corrplot(LbyE=LbyE,cor_type="heatmap",color=c("blue","white","red"))}
LbyE_corrplot <- function(LbyE,cor_type=NULL,color=NULL)
{
  if(is.null(cor_type)){
    cor_type = "heatmap"
  }
  if(is.null(color)){
    color=c("blue","white","red")
  }

  if(cor_type == "pie"){
    corrgram::corrgram(as.matrix(LbyE[,-1]), order=TRUE, lower.panel=panel.shade,
             pch = 19, upper.panel=panel.pie,
             col.regions = grDevices::colorRampPalette(c(color[1], color[2],high = color[3])));
  }else if(cor_type == "heatmap"){
    corr_data <- reshape2::melt(stats::cor(as.matrix(LbyE[,-1])))
    p<-ggplot(corr_data, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile() +
      geom_text(aes(label = round(value, 2)), color = "white") +
      scale_fill_gradient2(low = color[1], mid = color[2], high = color[3], midpoint = 0) +
      theme_minimal()+xlab(NULL)+ylab(NULL)
    print(p)
  }
}
