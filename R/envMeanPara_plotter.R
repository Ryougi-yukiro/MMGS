envMeanPara_plotter<-function(input,paras,method=NULL,linecolor=NULL,linewidth=NULL,
                              linetype=NULL,size=NULL,shape=NULL){
  input<-input[,c(1:5,match(paras, colnames(input)))]
  colnames(input)[ncol(input)]<-"Para"
  if(is.null(method)){
    method= "lm"
  }
  if(is.null(linecolor)){
    linecolor= "black"
  }
  if(is.null(linewidth)){
    linewidth= 1
  }
  if(is.null(linetype)){
    linetype= "dashed"
  }
  if(is.null(size)){
    size= 5
  }
  if(is.null(shape)){
    shape= 18
  }
  model <- lm(mean ~ Para, data = input)
  l <- list(r2 = format(summary(model)$r.squared, digits = 4))
  x1<-mean(input$Para);y1=min(input$mean);
  ggplot() +
    geom_point(data=input, aes(x = Para, y = mean ,color=env_code),size = size, shape = shape) +
    geom_smooth(data=input, aes(x = Para, y = mean ),method =method,
                se = FALSE, color=linecolor,linetype="dashed",linewidth=linewidth,formula = y ~ x)+
    theme_bw()+labs(x= paras, y= "mean")+
    geom_text(aes(x = x1, y = y1, label = (paste("R2",as.character(as.expression(l)),sep=' = '))))
}
