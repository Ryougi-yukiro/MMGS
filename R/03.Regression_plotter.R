#' @title Reg_plotter
#' @description
#'
#' Generate a regression plot between the mean and trait,
#' showing the relationship between these two variables.
#'
#' @param Reg Input data frame containing regression results.
#' param x Variable to be plotted on the x-axis.
#' param y Variable to be plotted on the y-axis.
#' param group Variable used for grouping the data.
#' @param size Size of points in the plot.
#' @param method Smoothing method for the regression line.
#' @param color Vector containing colors for regression line and points.
#' @param alpha Vector containing alpha values for regression line and points.
#'
#' @return A plot for the regression relationship between mean and trait.
#' @export
#'
#' @importFrom ggplot2 aes geom_smooth geom_point scale_fill_gradient2
#' labs theme_bw geom_text
#' @importFrom dplyr slice group_by  %>%
#' @examples
#' #' #Get Input
#' LbyE<-LbyE_calculate(trait,trait="FTgdd",env="env_code",line ="line_code")
#' env_trait<-env_trait_calculate(data=trait,trait ="FTgdd",env ="env_code")
#' Regression<-Reg(LbyE = LbyE, env_trait = env_trait)
#' Reg_plotter(Reg = Regression)
#'
Reg_plotter<-function(Reg=Reg,
                      size=NULL,
                      method=NULL,
                      color=NULL,
                      alpha=NULL){
  if(is.null(size)){
    size= 1
  }
  if(is.null(method)){
    method= "lm"
  }
  if(is.null(color)){
    color= c("gray90","gray50")
  }
  if(is.null(alpha)){
    alpha= c(0.1,0.5)
  }

  max_values <- Reg %>% group_by(Reg[,3]) %>% slice(which.min(Reg[,2]))
  max_values <-as.data.frame(max_values)
  ggplot()+geom_smooth(data=Reg, aes(x = Reg[,1], y = Reg[,2],group = Reg[,4]),
                       method = "lm", se = FALSE,color="gray")+
    geom_point(data=Reg, aes(x = Reg[,1],y = Reg[,2],group = Reg[,4]),
               color= "black")+
    geom_text(data=max_values,aes(x =max_values[,1], y=max_values[,2] ,label=max_values[,3]),
              hjust = -0.5, vjust = 0.5, color = "black", size = 3)+
    theme_bw()


}
