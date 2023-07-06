#' Mean_trait_plot
#'
#' @param Reg Reg result
#' @param MSE MSE reuslt
#' @param point_size point_size
#' @param point_shape point_shape
#' @param linewidth linewidth
#' @param linetype linetype
#'
#' @return
#' @export
#'
#' @examples Mean_trait_plot(Reg,MSE)
Mean_trait_plot<-function(Reg,MSE,point_size=NULL,
                          point_shape=NULL,linewidth=NULL,linetype=NULL){
  if(is.null(point_size)){
    point_size= 0.5
  }
  if(is.null(point_shape)){
    point_shape= 1
  }
  if(is.null(linewidth)){
    linewidth= 0.5
  }
  if(is.null(linetype)){
    linetype= 'dashed'
  }
  ggplot(Reg %>% left_join(MSE)) +
  geom_line(aes(meanY, Trait , group = line), color = "#80808023", show.legend = F)+
  geom_point(aes(meanY, Trait, group = line), color = "#80808023", show.legend = F) +
  geom_point(aes(meanY, Trait, group = line, color = env_code),size = point_size, shape = point_shape) +
  geom_smooth(aes(meanY, Trait), color = 'grey', linewidth = linewidth,
              method = 'lm', se = F, linetype = linetype) +
  geom_ribbon(aes(meanY, ymin = q25, ymax = q75), fill = "#EE82EE37") +
  geom_point(aes(meanY, meanY))+theme_bw()
}
