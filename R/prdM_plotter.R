#' prdM_plotter
#'
#' @param prdM prdM
#' @param data envMeanPara
#' @param trait trait
#' @param Para_Name env para name
#' @param size point size
#' @param shape point shape
#' @param alpha point alpha
#' @param method regression method
#' @param linewidth line width
#' @param linetype line type
#' @param linecolor line color
#'
#' @param x_break adjust x label text
#' @param y_break adjust y label text
#'
#' @return combined prdM plot
#' @export
#'
#' @examples \dontrun{prdM_plotter(prdM=prdM,data=envMeanPara,trait="FTgdd",Para_Name="PTT")}
  prdM_plotter<-function(prdM,data,trait,Para_Name,size=NULL,shape=NULL,alpha=NULL,method=NULL,
                         linewidth=NULL,linetype=NULL,linecolor=NULL,x_break=NULL,y_break=NULL){
    if(is.null(size)){
      size= 1
    }
    if(is.null(shape)){
      shape= 20
    }
    if(is.null(method)){
      method= "lm"
    }
    if(is.null(linewidth)){
      linewidth= 1
    }
    if(is.null(alpha)){
      alpha= 1
    }
    if(is.null(linetype)){
      linetype= "dashed"
    }
    if(is.null(linecolor)){
      linecolor= "black"
    }
    if(is.null(x_break)){
      x_break= 100
    }
    if(is.null(y_break)){
      y_break= 20
    }
    prdM <- prdM[!is.na(prdM$Obs_trait),];
    prd_env <- as.vector(unique(prdM$env_code));
    env_rs <- matrix(ncol = 3, nrow = length(prd_env));
    for (e_i in 1:length(prd_env)) {
      env_obs_prd <- subset(prdM, prdM$env_code == prd_env[e_i]);
      if (nrow(env_obs_prd) > 0) {
        env_rs[e_i,] <- c(
          sprintf( "%.2f", cor(as.numeric(env_obs_prd[,4]), as.numeric(env_obs_prd[,6]),
                               use = "complete.obs")),
          sprintf( "%.2f", cor(as.numeric(env_obs_prd[,5]), as.numeric(env_obs_prd[,6]),
                               use = "complete.obs")),
          sprintf( "%.2f", cor(as.numeric(env_obs_prd[,7]), as.numeric(env_obs_prd[,6]),
                               use = "complete.obs")));
      }
    }
    colnames(env_rs)<-c("Prd_trait_mean ","Prd_trait_kPara","Line_mean")
    rownames(env_rs)<-prd_env
    #return(env_rs)
    prdM$Obs_trait<-as.numeric(prdM$Obs_trait)
    prdM$Prd_trait_mean<-as.numeric(prdM$Prd_trait_mean)
    prdM$Prd_trait_kPara<-as.numeric(prdM$Prd_trait_kPara)
    prdM$Line_mean<-as.numeric(prdM$Line_mean)

    x1<-(max(prdM$Obs_trait)+min(prdM$Obs_trait))/2;
    x2<-(max(data[[Para_Name]])+min(data[[Para_Name]]))/2;
    x3<-min(prdM$Obs_trait)+(max(prdM$Obs_trait)+min(prdM$Obs_trait))/x_break;
    y1=max(prdM$Prd_trait_mean);
    y1in <- (max(prdM$Prd_trait_mean)-min(prdM$Prd_trait_mean))/y_break;
    y2=max(prdM$Prd_trait_kPara);
    y2in <- (max(prdM$Prd_trait_kPara)-min(prdM$Prd_trait_kPara))/y_break;
    y3=max(prdM$Line_mean);
    y3in <- (max(prdM$Line_mean)-min(prdM$Line_mean))/y_break;
    y4=max(data$mean);


    colnames(data)[colnames(data) == Para_Name]<-"Paras"

    model <- lm(Obs_trait ~ Prd_trait_mean, data = prdM)
    l1 <- list(r2 = format(summary(model)$r.squared, digits = 4))

    d1<-data.frame(env_code = character(), R2 = numeric(), stringsAsFactors = FALSE)
    for (i in unique(prdM$env_code)){
      data1<- subset(prdM, env_code == i)
      r2<-summary(lm(formula = Obs_trait ~ Prd_trait_mean, data = data1))$r.squared
      d1 <- rbind(d1, data.frame(env_code = i, R2 = r2, stringsAsFactors = FALSE))
    }

    model <- lm(Obs_trait ~ Prd_trait_kPara, data = prdM)
    l2 <- list(r2 = format(summary(model)$r.squared, digits = 4))

    d2<-data.frame(env_code = character(), R2 = numeric(), stringsAsFactors = FALSE)
    for (i in unique(prdM$env_code)){
      data2<- subset(prdM, env_code == i)
      r2<-summary(lm(formula = Obs_trait ~ Prd_trait_kPara, data = data2))$r.squared
      d2 <- rbind(d2, data.frame(env_code = i, R2 = r2, stringsAsFactors = FALSE))
    }

    model <- lm(Obs_trait ~ Line_mean, data = prdM)
    l3 <- list(r2 = format(summary(model)$r.squared, digits = 4))
    d3<-data.frame(env_code = character(), R2 = numeric(), stringsAsFactors = FALSE)
    for (i in unique(prdM$env_code)){
      data3<- subset(prdM, env_code == i)
      r2<-summary(lm(formula = Obs_trait ~ Line_mean, data = data3))$r.squared
      d3<- rbind(d3, data.frame(env_code = i, R2 = r2, stringsAsFactors = FALSE))
    }

    model <- lm(Paras ~ mean, data = data)
    l4 <- list(r2 = format(summary(model)$r.squared, digits = 4))

    p1<-ggplot()+geom_point(data=prdM,aes(x=Obs_trait,y=Prd_trait_mean,color=env_code),size=size,
                            shape=shape,alpha=alpha)+
      geom_smooth(data=prdM,aes(x=Obs_trait,y=Prd_trait_mean),se=F,linewidth=linewidth,method=method,
                  linetype=linetype,color=linecolor,formula = y ~ x)+theme_bw()+
      labs(x=paste("Observed",trait,sep=" "),
           y=paste("Predicted by","envMean",sep=" "))+
      geom_text(aes(x = x1, y = y1,
                    label = (paste("R2",as.character(as.expression(l1)),sep=' = '))))
    for(i in 1:nrow(d1)){
      p1<-p1+annotate('text',x = x3,y = y1-y1in*(i-1),
                      label = (paste("R2 =", round(d1[i,2], 2)," (",d1[i,1],")")),
                      size=3,vjust = 0.5,hjust = 0.5)
    }

    p2<-ggplot()+geom_point(data=prdM,aes(x=Obs_trait,y=Prd_trait_kPara,color=env_code),size=size,
                            shape=shape,alpha=alpha)+
      geom_smooth(data=prdM,aes(x=Obs_trait,y=Prd_trait_kPara),se=F,linewidth=linewidth,method=method,
                  linetype=linetype,color=linecolor,formula = y ~ x)+theme_bw()+
      labs(x=paste("Observed",trait,sep=" "),
           y=paste("Predicted by",Para_Name,sep=" "))+
      geom_text(aes(x = x1, y = y2,
                    label = (paste("R2",as.character(as.expression(l2)),sep=' = '))))
    for(i in 1:nrow(d2)){
      p2<-p2+annotate('text',x = x3,y = y2-y2in*(i-1),
                      label = (paste("R2 =", round(d2[i,2], 2)," (",d2[i,1],")")),
                      size=3,vjust = 0.5,hjust = 0.5)
    }

    p3<-ggplot()+geom_point(data=prdM,aes(x=Obs_trait,y=Line_mean,color=env_code),size=size,
                            shape=shape,alpha=alpha)+
      geom_smooth(data=prdM,aes(x=Obs_trait,y=Line_mean),se=F,linewidth=linewidth,method=method,
                  linetype=linetype,color=linecolor,formula = y ~ x)+theme_bw()+
      labs(x=paste("Observed",trait,sep=" "),
           y=paste("Predicted by","BLUE",sep=" "))+
      geom_text(aes(x = x1, y = y3,
                    label = (paste("R2",as.character(as.expression(l3)),sep=' = '))))
    for(i in 1:nrow(d3)){
      p3<-p3+annotate('text',x = x3,y = y3-y3in*(i-1),
                      label = (paste("R2 =", round(d3[i,2], 2)," (",d3[i,1],")")),
                      size=3,vjust = 0.5,hjust = 0.5)
    }
    p4<-ggplot()+geom_point(data=data,aes(x=Paras,y=mean,color=env_code),size=size,
                            shape=shape,alpha=alpha)+
      geom_smooth(data=data,aes(x=Paras,y=mean),se=F,linewidth=linewidth,method=method,
                  linetype=linetype,color=linecolor,formula = y ~ x)+theme_bw()+
      labs(x=Para_Name,y="Observed Mean")+
      geom_text(aes(x = x2, y = y4,
                    label = (paste("R2",as.character(as.expression(l4)),sep=' = '))))

    lemon::grid_arrange_shared_legend(p1, p2, p3, p4, nrow = 2, ncol = 2,position='right')
  }
