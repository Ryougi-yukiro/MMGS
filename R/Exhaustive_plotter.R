#' Exhaustive_plotter
#'
#' @param data pop_cor file
#' @param dap_x depend on your file
#' @param dap_y depend on your file
#' @param p defualt p=1
#'
#' @return combined pop_cor plot
#' @export
#'
#' @examples \dontrun{Exhaustive_plotter(input=pop_cor,paras="PTS")}
Exhaustive_plotter<-function(data,dap_x, dap_y,p){
pop_cors<-as.data.frame(data)
ps<-list()
for (k in 1:(2*length(Paras))) {
  pop_cor <- pop_cors[,c(1:4, k + 4)];
  title_text <- colnames(pop_cor)[ncol(pop_cor)]
  colnames(pop_cor)[ncol(pop_cor)] <- "Cor"
  Max_cor<- pop_cor[which.max(pop_cor$Cor), ]
  p<-ggplot(pop_cor, aes(x = Day_x, y =  Day_y, fill = Cor))+geom_tile()+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
    geom_text(data = Max_cor, aes(label =  paste("x:", `Day_x`, "\ny:", `Day_y`,"\ncor:",`Cor`)),
              color = "black", size = 4)+
    geom_label(data = Max_cor, aes(label = paste("x:", Day_x, "\ny:", Day_y, "\ncor:", Cor)),
               color = "black", size = 4, fill = "yellow", alpha = 0.8, show.legend = FALSE)+
    theme(axis.ticks = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),)+
    theme_bw()+scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(expand=c(0,0))+labs(title = title_text,vjust=2)

  ps[[k]] <- p
}
d <- character(length(ps))
for (i in 1:length(ps)) {d[i] <- paste0("ps[[", i,"]]")}
a<-paste(d, collapse = ", ")
code <- paste("lemon::grid_arrange_shared_legend(", paste(a, collapse = ", "),
              ", nrow = 2, ncol = 5, position = 'right')")
eval(parse(text = code))
}
