
#' Predict
#'
#' @param pheno pheno info, from LbyE
#' @param geno geno info
#' @param env env info
#' @param para envMeanPara calculated by envMeanPara
#' @param Para_Name Your env para Name
#' @param depend Maker or Norm
#' @param fold Number of folds like rrblup
#' @param reshuffle Repeat the calculation for how many cycles
#' @param methods RM.G RM.GE RM.E
#' @param ms2 remove the line for which the number of missing environment > ms2
#' @param line_missing remove the environment for which the number of missing line > line_missing
#'
#' @return prediction result
#' @export
#'
#' @examples out<-Predict(pheno=pheno,geno=geno,env=env_info,para=envMeanPara,Para_Name="PTT",
#'                        depend="maker",fold=2,reshuffle=5,methods="RM.G",
#'                        ms1=ms1,ms2=ms2)
Predict<-function(pheno,geno,env,para,Para_Name,depend=NULL,fold=NULL,reshuffle=NULL,
                  methods=NULL,ms1=NULL,ms2=NULL){
  if(is.null(depend)){
    depend="norm"
  }
  if(is.null(methods)){
    methods="RM.GE"
  }
  if(is.null(fold)){
    fold= 5
  }
  if(is.null(reshuffle)){
    reshuffle= 5
  }
  if(is.null(ms2)){
    ms2= c(ncol(pheno)-1)-3;
  }
  if(is.null(ms1)){
    ms1= nrow(pheno)*0.5;
  }
  enp=which(colnames(para) == Para_Name);
  para.name=colnames(pheno)[-1];
  #m=pheno[,-1];
  #m=as.numeric(as.character(unlist(m)));
  #m <- matrix(data=m, ncol=dim(pheno)[2]-1, nrow=dim(pheno)[1]);
  #colnames(m)=colnames(pheno)[-1];
  #pheno_=data.frame(line_code=pheno$line_code,m);
  #colnames(pheno_)=c("line_code",para.name);
  #pheno=pheno_;

  para=para[para$env_code%in%para.name,];
  env=env[env$env_code%in%para.name,];
  order=match(para$env_code,env$env_code);
  env=env[order,];
  rm.env=colnames(pheno)[colSums(is.na(pheno))>=ms1];

  pheno=pheno[,which(!(colnames(pheno)%in%rm.env))];
  env=env[which(!(env$env_code%in%rm.env)),];
  para=para[which(!(para$env_code%in%rm.env)),];

  pheno.1=pheno
  keep=dim(pheno.1)[2]-rowSums(is.na(pheno.1))
  pheno=pheno.1[which(keep>=ms2),];
  n.line=dim(pheno)[1];
  n.para=dim(pheno)[2]-1;

  coloo=heat.colors(n.para);
  if(depend == "norm"){
  if(methods=="RM.E"){
    pheno.hat=matrix(999,n.line,dim(para)[1]);
    cor.whole=numeric();
    for(k in 1:n.para){
      for(j in 1:n.line){
        x1=para[,enp][-k];
        y1=as.vector(t(pheno[j,-c(1,1+k)]));
        coe=lm(y~x,data=data.frame(x=x1,y=y1));
        y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*para[,enp][k];
        pheno.hat[j,k]=y.hat;
      }
      cor.each=cor(as.numeric(pheno.hat[,k]),as.numeric(pheno[,k+1]), use = "complete.obs");
      cor.whole=c(as.numeric(cor.whole),as.numeric(cor.each));
    }
    observe=as.vector(as.matrix(pheno[,-1]));
    predict=as.vector(pheno.hat);

    r_within=cor.whole;y1
    names(r_within)=colnames(pheno)[-1];
    r_within=data.frame(cor_within=r_within,para=colnames(pheno)[-1]);
    r_across=cor(as.numeric(observe),as.numeric(predict),use = "complete.obs");
    outforfigure=data.frame(obs=observe,pre=predict,
                            col=rep(coloo,times=rep(n.line,n.para)),
                            para=rep(colnames(pheno)[-1],times=rep(n.line,n.para)))
    #return(outforfigure)
    return(list(outforfigure,r_within,r_across));
    #return(r_within);
  }
  if(methods=="RM.G"){
    intercept=numeric();
    slope=numeric();
    for(j in 1:n.line)
    {
      x1=para[,enp];
      y1=as.vector(t(pheno[j,-c(1)]));

      coe=lm(y~x,data=data.frame(x=x1,y=y1));
      inter=summary(coe)$coefficients[1,1]
      slop=summary(coe)$coefficients[2,1];
      intercept=c(intercept,inter);
      slope=c(slope,slop);
    }
    genotype.match=match(pheno[,1],geno[,1])
    genotype_1=geno[genotype.match,];
    genotype.impute=A.mat(genotype_1[,-1],max.missing=0.5,impute.method="mean",return.imputed=T);
    SFI=cbind(genotype_1[,1],genotype.impute$imputed);
    genotype=matrix(suppressWarnings(as.numeric(SFI)),nrow=nrow(SFI))
    Marker=genotype[,-1];
    intercept.hat=numeric();
    slope.hat=numeric();
    cor.within=matrix(999,reshuffle,n.para);
    cor.all=numeric();
    for(i in 1:reshuffle){
      cross = sample(rep(1:fold,each=ceiling(n.line/fold)),n.line)
      yhat.whole.cross=numeric();yobs.whole.cross=numeric();
      for(f in 1:fold){
        id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
        y0=intercept;
        ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
        e=as.matrix(ans$u)
        G.pred=Marker[id.V,]
        y_pred=as.matrix(G.pred) %*% e
        GEBV.inter=c(y_pred[,1])+c(ans$beta);
        y0=slope;
        ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
        e=as.matrix(ans$u)
        G.pred=Marker[id.V,]
        y_pred=as.matrix(G.pred) %*% e
        GEBV.slope=c(y_pred[,1])+c(ans$beta);
        yhat.para=matrix(999,length(id.V),n.para);yobs.para=matrix(999,length(id.V),n.para)
        for(j in 1:n.para)
        {
          yhat=GEBV.inter+GEBV.slope*para[j,enp];
          yobs=pheno[id.V,j+1];
          yhat.para[,j]=yhat;yobs.para[,j]=yobs;
        }
        yhat.whole.cross=rbind(yhat.whole.cross,yhat.para);
        yobs.whole.cross=rbind(yobs.whole.cross,yobs.para);
      }
      for(j in 1:n.para)
      {cor.within[i,j]=cor(yhat.whole.cross[,j],yobs.whole.cross[,j],use = "complete.obs");}

      cor.all=c(cor.all,cor(as.vector(yhat.whole.cross),
                            as.vector(yobs.whole.cross),use = "complete.obs"));
    }
    #Correlation within environment 50 times
    r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1];
    r_within=data.frame(cor_within=r_within,para=colnames(pheno)[-1]);
    #Correlation across environment 50 times
    r_across=mean(cor.all);
    #Observation and prediction last time
    outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                            pre=as.vector(yhat.whole.cross),
                            col=rep(coloo,times=rep(n.line,n.para)),
                            para=rep(colnames(pheno)[-1],times=rep(n.line,n.para)))
    colnames(cor.within)=colnames(pheno)[-1];
    r_within=cor.within;
    r_across=cor.all;
    return(list(outforfigure,r_within,r_across));
      }
  if(methods=="RM.GE"){
    genotype.match=match(pheno[,1],geno[,1])
    genotype_1=geno[genotype.match,];
    genotype.impute=A.mat(genotype_1[,-1],max.missing=0.5,impute.method="mean",return.imputed=T);
    SFI=cbind(genotype_1[,1],genotype.impute$imputed);
    genotype=matrix(suppressWarnings(as.numeric(SFI)),nrow=nrow(SFI))
    Marker=genotype[,-1];

    cor.within=matrix(999,reshuffle,n.para);cor.all=numeric();
    for(i in 1:reshuffle)
    {
      obs_matrix=matrix(999,n.line,n.para);pre_matrix=matrix(999,n.line,n.para);
      for(k in 1:n.para)
      {
        intercept=numeric();
        slope=numeric();
        for(j in 1:n.line)
        {
          x1=para[-k,enp];
          y1=as.vector(t(pheno[j,-c(1,1+k)]));

          coe=lm(y~x,data=data.frame(x=x1,y=y1));
          inter=summary(coe)$coefficients[1,1]
          slop=summary(coe)$coefficients[2,1];
          intercept=c(intercept,inter);
          slope=c(slope,slop);
        }

        cross = sample(rep(1:fold,each=ceiling(n.line/fold)),n.line)
        yhat.whole=numeric();
        yobs.whole=numeric();

        for(f in 1:fold)
        {
          id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
          y0=intercept;
          ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
          e=as.matrix(ans$u)
          G.pred=Marker[id.V,]
          y_pred=as.matrix(G.pred) %*% e
          GEBV.inter=c(y_pred[,1])+c(ans$beta);
          y0=slope;
          ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
          e=as.matrix(ans$u)
          G.pred=Marker[id.V,]
          y_pred=as.matrix(G.pred) %*% e
          GEBV.slope=c(y_pred[,1])+c(ans$beta);
          ###All the predicted slope and intercept
          yhat=GEBV.inter+GEBV.slope*para[k,enp];
          yobs=pheno[id.V,k+1];

          yhat.whole=c(yhat.whole,yhat);
          yobs.whole=c(yobs.whole,yobs);
        }
        cor.within[i,k]=cor(yhat.whole,yobs.whole,use = "complete.obs");
        obs_matrix[,k]=yobs.whole;
        pre_matrix[,k]=yhat.whole;
      }
      cor.shuffle=cor(as.vector(obs_matrix),as.vector(pre_matrix),use = "complete.obs")
      cor.all=c(cor.all,cor.shuffle);
    }
    yhat.whole.cross=pre_matrix;
    yobs.whole.cross=obs_matrix;
    r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1];
    r_across=mean(cor.all);
    outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                            pre=as.vector(yhat.whole.cross),
                            col=rep(coloo,times=rep(n.line,n.para)),
                            para=rep(colnames(pheno)[-1],times=rep(n.line,n.para)))
    colnames(cor.within)=colnames(pheno)[-1];
    r_within=cor.within;
    r_across=cor.all;
    return(list(outforfigure,r_within,r_across));
  }
  }else if(depend =="maker"){
    genotype.match=match(pheno[,1],geno[,1])
    genotype_1=geno[genotype.match,];
    genotype.impute=A.mat(genotype_1[,-1],max.missing=0.5,impute.method="mean",return.imputed=T);
    SFI=cbind(genotype_1[,1],genotype.impute$imputed);
    genotype=matrix(suppressWarnings(as.numeric(SFI)),nrow=nrow(SFI))
    Marker=genotype[,-1];
    n.marker=dim(Marker)[2];
    if(methods=="RM.E"){
      effect=matrix(999,n.marker,n.para);
      intercept=matrix(999,1,n.para)
      for(i in 1:n.para)
      {
        fit=mixed.solve(pheno[,1+i],Z=Marker)
        effect[,i]=fit$u
        intercept[,i]=fit$beta
      }
      pheno.hat=matrix(999,n.line,dim(para)[1]);
      cor.whole=numeric();
      for(k in 1:n.para)
      {
        effect.hat=numeric();
        for(j in 1:n.marker)
        {
          x1=para[,enp][-k];
          y1=effect[j,-k];

          coe=lm(y~x,data=data.frame(x=x1,y=y1));
          y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*para[,enp][k];
          effect.hat=c(effect.hat,y.hat);
        }
        reg.intercept=data.frame(x=as.numeric(para[,enp][-k]),y=as.vector(intercept[-k]));
        coe.intercept=lm(y~x,data=reg.intercept)
        y.intercept=summary(coe.intercept)$coefficients[1,1]+summary(coe.intercept)$coefficients[2,1]*as.numeric(para[,enp][k]);

        pheno.hat[,k]=y.intercept+as.matrix(Marker)%*%as.matrix(effect.hat);
        cor.each=cor(pheno.hat[,k],pheno[,k+1], use = "complete.obs");
        cor.whole=c(cor.whole,cor.each);
      }
      observe=as.vector(as.matrix(pheno[,-1]));
      predict=as.vector(pheno.hat);
      r_within=cor.whole;names(r_within)=colnames(pheno)[-1];
      r_within=data.frame(cor_within=r_within,para=colnames(pheno)[-1]);
      r_across=cor(observe,predict,use = "complete.obs");
      outforfigure=data.frame(obs=observe,pre=predict,col=rep(coloo,times=rep(n.line,n.para)),
                              para=rep(colnames(pheno)[-1],times=rep(n.line,n.para)))
      return(list(outforfigure,r_within,r_across));
    }
    if(methods=="RM.G"){
      cor.within=matrix(999,reshuffle,n.para);
      cor.all=numeric();
      for(i in 1:reshuffle)
      {
        cross = sample(rep(1:fold,each=ceiling(n.line/fold)),n.line)
        yhat.whole.cross=numeric();
        yobs.whole.cross=numeric();
        for(f in 1:fold)
        {
          id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
          ##marker effect###
          effect=matrix(999,n.marker,n.para);
          intercept=matrix(999,1,n.para)

          for(k in 1:n.para)
          {
            fit=mixed.solve(pheno[id.T,1+k],Z=Marker[id.T,])
            effect[,k]=fit$u
            intercept[,k]=fit$beta
          }

          ##Slope###
          effect.hat=numeric();
          for(j in 1:n.marker)
          {
            x1=para[,enp];
            y1=effect[j,];
            coe=lm(y~x,data=data.frame(x=x1,y=y1));
            y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*para[,enp];
            effect.hat=rbind(effect.hat,y.hat);
          }

          ##Environment intercept####
          reg.intercept=data.frame(x=as.numeric(para[,enp]),y=as.vector(intercept));
          coe.intercept=lm(y~x,data=reg.intercept)
          y.intercept=summary(coe.intercept)$coefficients[1,1]+
            summary(coe.intercept)$coefficients[2,1]*as.numeric(para[,enp]);

          ###All the predicted slope and intercept
          yhat.para=matrix(999,length(id.V),n.para);
          yobs.para=matrix(999,length(id.V),n.para)
          for(j in 1:n.para)
          {
            yhat=y.intercept[j]+(as.matrix(Marker[id.V,])%*%as.matrix(effect.hat))[,j];
            yobs=pheno[id.V,j+1];
            yhat.para[,j]=yhat;
            yobs.para[,j]=yobs;
          }
          yhat.whole.cross=rbind(yhat.whole.cross,yhat.para);
          yobs.whole.cross=rbind(yobs.whole.cross,yobs.para);
        }

        for(j in 1:n.para)
        {cor.within[i,j]=cor(yhat.whole.cross[,j],yobs.whole.cross[,j],use = "complete.obs");}
        cor.all=c(cor.all,cor(as.vector(yhat.whole.cross),
                              as.vector(yobs.whole.cross),use = "complete.obs"));
      }
      #Correlation within environment 50 times
      r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1];
      #Correlation across environment 50 times
      r_across=mean(cor.all);
      #Observation and prediction last time
      outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                              pre=as.vector(yhat.whole.cross),
                              col=rep(coloo,times=rep(n.line,n.para)),
                              para=rep(colnames(pheno)[-1],times=rep(n.line,n.para)))
      colnames(cor.within)=colnames(pheno)[-1];
      r_within=cor.within;
      r_across=cor.all;
      return(list(outforfigure,r_within,r_across))
    }
    if(methods=="RM.GE"){
      cor.within=matrix(999,reshuffle,n.para);cor.all=numeric();
      for(i in 1:reshuffle)
      {
        cross = sample(rep(1:fold,each=ceiling(n.line/fold)),n.line)
        obs_matrix=numeric();pre_matrix=numeric();
        for(f in 1:fold)
        {
          id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
          effect=matrix(999,n.marker,n.para);
          intercept=matrix(999,1,n.para)

          for(k in 1:n.para)
          {
            fit=mixed.solve(pheno[id.T,1+k],Z=Marker[id.T,])
            effect[,k]=fit$u
            intercept[,k]=fit$beta
          }

          ##predict marker effect###
          obs.para=numeric();pre.para=numeric();
          for(kk in 1:n.para)
          {
            effect.hat=numeric();
            for(j in 1:n.marker)
            {
              x1=para[-kk,enp];
              y1=effect[j,-kk];
              coe=lm(y~x,data=data.frame(x=x1,y=y1));
              y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*para[kk,enp];
              effect.hat=rbind(effect.hat,y.hat);
            }
            ##Environment intercept####
            reg.intercept=data.frame(x=as.numeric(para[-kk,enp]),y=as.vector(intercept[-kk]));
            coe.intercept=lm(y~x,data=reg.intercept)
            y.intercept=summary(coe.intercept)$coefficients[1,1]+summary(coe.intercept)$coefficients[2,1]*as.numeric(para[kk,enp]);

            ###All the predicted slope and intercept
            yhat=y.intercept+(as.matrix(Marker[id.V,])%*%as.matrix(effect.hat));
            yobs=pheno[id.V,kk+1];
            obs.para=cbind(obs.para,yobs);
            pre.para=cbind(pre.para,yhat);
          }
          obs_matrix=rbind(obs_matrix,obs.para);pre_matrix=rbind(pre_matrix,pre.para);

        }
        cor.shuffle=cor(as.vector(obs_matrix),as.vector(pre_matrix),use = "complete.obs")
        cor.all=c(cor.all,cor.shuffle);
        for(j in 1:n.para)
        {cor.within[i,j]=cor(obs_matrix[,j],pre_matrix[,j],use = "complete.obs");}
      }
      yhat.whole.cross=pre_matrix;
      yobs.whole.cross=obs_matrix;
      #Correlation within environment 50 times
      r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1];
      #Correlation across environment 50 times
      r_across=mean(cor.all);
      #Observation and prediction last time
      outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                              pre=as.vector(yhat.whole.cross),
                              col=rep(coloo,times=rep(n.line,n.para)),para=rep(colnames(pheno)[-1],times=rep(n.line,n.para)))
      colnames(cor.within)=colnames(pheno)[-1];
      r_within=cor.within;
      r_across=cor.all;
      return(list(outforfigure,r_within,r_across));
    }
  }
}
