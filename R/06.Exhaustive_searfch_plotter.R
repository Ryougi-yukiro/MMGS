#' Eh_plot
#'
#' Generate a combined plot for pairwise correlations using pop_cor data.
#'
#' @param Correlation A pop_cor data frame containing correlation values for different pairs.
#' @param dap_x Specifies the column for the x-axis variable.
#' @param dap_y Specifies the column for the y-axis variable.
#' @param p A parameter with a default value of 1.
#' @param Paras Vector containing the names of environmental parameters to be considered.
#' @param Cor Assigning variable column names.Default is Cor.
#' @return A combined pop_cor plot.
#' @export
#'
#' @importFrom ggplot2 theme geom_label element_blank theme_bw
#' scale_y_continuous scale_x_continuous labs
#' @importFrom lemon grid_arrange_shared_legend
#' @examples
#' env_trait<-env_trait_calculate(data=trait,trait="FTgdd",env="env_code")
#' #Run
#' \donttest{pop_cor<-Exhaustive_search(env_trait,PTT_PTR,80)
#' dap1<-pop_cor[,2][which(abs(pop_cor[, -c(1:4)]) ==
#' max(abs(pop_cor[, -c(1:4)])), arr.ind = TRUE)][1]
#' dap2<-pop_cor[,3][which(abs(pop_cor[, -c(1:4)]) ==
#' max(abs(pop_cor[, -c(1:4)])), arr.ind = TRUE)][1]
#' Eh_plot(pop_cor,dap1,dap2)}
Eh_plot <- function(Correlation,
                    dap_x,dap_y,
                    p = NULL,
                    Paras=NULL,Cor = NULL) {
  if(is.null(Cor)){Cor=Cor}
  if(is.null(p)){p=1}
  if(is.null(Paras)){
    Paras= c("DL", "GDD", "PTT", "PTR", "PTS")
  }
  pop_cors <- as.data.frame(Correlation)
  ps <- list()
  # Loop through the columns of the pop_cor data
  for (k in 1:(2 * length(Paras))) {
    Max_cor<-data.frame()
    df1<-data.frame()
    # Extract relevant columns for the correlation plot
    df1<-data.frame(
      pop_code = pop_cors[,1],
      Day_x = pop_cors[,2],
      Day_y = pop_cors[,3],
      window = pop_cors[,4]
    )
    df1<-as.data.frame(df1)
    df1$Cor <-pop_cors[,k+4]
    #colnames(df1)[5] <- colnames(pop_cors)[k+4]
    #print(head(df1))
    #pop_cor <- pop_cors[, c(1:4, k + 4)]
    #title_text <- colnames(df1)[ncol(df1)]
    #colnames(df1)[ncol(df1)] <- "Cor"
    Max_cor <- df1[which.max(df1$Cor), ]
    # Create ggplot for the correlation plot
    p <- ggplot()+
      geom_tile(data=df1,
                aes(x = df1[,2], y = df1[,3], fill = Cor)) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0, limits = c(-1, 1)) +
      #geom_point(data=Max_cor,
      #           aes(x=data[,2],y=data[,3],color="purple",size=3))+
      annotate("text", x = dap_x-15, y = 10,
               label = paste("x:", Max_cor[,2], "\ny:",
                             Max_cor[,3],"\nCor:",Max_cor[,5]),
               hjust = 1, vjust = 0)+
      #geom_label(data = Max_cor, aes(label = paste("x:", Max_cor[,2], "\ny:",
      #                                             Max_cor[,3], "\ncor:",
      #                                              Max_cor[,5]
      #           color = "black", size = 4, fill = "pink",
      #            alpha = 0.8, show.legend = FALSE) +
      theme(axis.ticks = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),) +
      theme_bw() +
      xlab("Dap_x")+ylab("Dap_y")+
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(title =  colnames(pop_cors)[k+4], vjust = 2,
           fill =  "Cor")

    ps[[k]] <- p
  }

  # Combine plots into a grid with shared legend
  d <- character(length(ps))
  for (i in 1:length(ps)) {
    d[i] <- paste0("ps[[", i, "]]")
  }
  a <- paste(d, collapse = ", ")
  code <- paste("lemon::grid_arrange_shared_legend(", paste(a, collapse = ", "),
                ", nrow = 2, ncol = 5, position = 'right')")
  eval(parse(text = code))
}
