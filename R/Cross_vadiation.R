#' Cross Vadidation and Genomic Prediction Cross environments
#' @description
#' The GP_CV function carries out cross-validation using genotypic and phenotypic data from a reference population,
#' with result for Genomic Prediction for different environments,  genomic breeding value estimation for each env and
#' cross R2 for all envs.
#'
#' @param geno Matrix (n x m) of genotypes for the training population: n lines with m markers.
#'             Genotypes should be coded -1, 0, 1. Missing data are not allowed.
#' @param pheno Vector (n x j) of "phenotypes", i.e. observations or pre-processed, corrected values.
#' @param env Data.frame (j x l) of "environmental information", i.e.
#' @param para Data.frame returned from envMeanPara function.
#' @param Para_Name The most relevant environmental silvers to the subject's phenotype,
#'                  obtained after stepwise correlation calculations, are referred to for more details:
#' @param depend The options for genomic breeding within different environment.The available options are:
#'               1.Norm:
#'               2.Marker
#' @param fold 	Number of folds for the cross-validation.
#'              Smallest value recommended is fold = 2.
#' @param reshuffle Number of independent replicates for the cross-validation.
#'                  Smallest value recommended is reshuffle = 5.
#' @param methods RM.G RM.GE RM.E
#' @param model The options for genomic breeding cross vadiation methods. The available options are:
#' 1.GBLUP: performs G-BLUP using a marker-based relationship matrix, implemented through BGLR R-library.
#'              Equivalent to ridge regression (RR-BLUP) of marker effects.
#' 2.RR: ridge regression, using package glmnet.
#'              In theory, strictly equivalent to gblup.
#' 3.LASSO: Least Absolute Shrinkage and Selection Operator is another penalized regression methods
#'              which yield more shrinked estimates than RR.Run by glmnet library.
#'              4.EN: Elastic Net (Zou and Hastie, 2005),
#'              which is a weighted combination of RR and LASSO, using glmnet library
#'              5.rrBLUP:
#'              Several Bayesian methods, using the BGLR library:
#'              1.BRR: Bayesian ridge regression: same as rr-blup, but bayesian resolution.
#'              Induces homogeneous shrinkage of all markers effects towards zero with
#'              Gaussian distribution (de los Campos et al, 2013)
#'              2.BL: Bayesian LASSO: uses an exponential prior on marker variances priors,
#'              leading to double exponential distribution of marker effects (Park & Casella 2008)
#'              3.BA: uses a scaled-t prior distribution of marker effects. (Meuwissen et al 2001).
#'              4.BB: Bayes B, uses a mixture of distribution with a point mass at zero
#'              and with a slab of non-zero marker effects with a scaled-t distribution (Habier et al 2011).
#'              5.BC: Bayes C same as Bayes B with a slab with Gaussian distribution.
#'              A more detailed description of these methods can be found in Perez & de los Campos 2014 (http://genomics.cimmyt.org/BGLR-extdoc.pdf).
#'              Four semi-parametric methods:
#'              1.RKHS: reproductive kernel Hilbert space and multiple kernel MRKHS, using BGLR (Gianola and van Kaam 2008).
#'              Based on genetic distance and a kernel function to regulate the distribution of marker effects.
#'              This methods is claimed to be effective for detecting non additive effects.
#'              2.RF: Random forest regression, using randomForest library (Breiman, 2001, Breiman and Cutler 2013).
#'              This methods uses regression models on tree nodes which are rooted in bootstrapping data.
#'              Supposed to be able to capture interactions between markers
#'              3.SVM: support vector machine, run by e1071 library. For details,
#'              see Chang, Chih-Chung and Lin, Chih-Jen: LIBSVM: a library for Support Vector Machines http://www.csie.ntu.edu.tw/~cjlin/libsvm
#'              4.LightGBM: Light Gradient Boosting Machine, run by LightGBM library developed by Microsoft.
#'              For details, see Jun Yan, Yuetong Xu, at al: LightGBM: accelerated genomically designed crop breeding through ensemble learning
#'              and Mircosoft Mannual https://lightgbm.readthedocs.io/en/latest/R/index.html
#' @param ms1 remove the line for which the number of missing environment > ms1
#' @param ms2 remove the line for which the number of missing environment > ms2
#' @param alpha used for EN predict model,the value range from 0 to 1.
#' @param nIter Number of iterations for the Bayesian model.
#' @param burnIn Number of burn-in for the Bayesian model.
#' @param thin Number of thins for the Bayesian model.
#' @param SVM_cost cost of constraints violation (default: 10)—it is the 'C'-constant of the regularization term in the Lagrange formulation
#' @param gamma parameter needed for all kernels except linear (default: 0.001)
#' @param GBM_params  a list of parameters. See The "Dataset Parameters"{https://lightgbm.readthedocs.io/en/latest/Parameters.html#dataset-parameters} section of the documentation
#' @param GBM_rounds  number of training rounds Based on LightGBM model
#'
#' @return The class GP_CV returns a list containing:
#'         The most relevant environmental silvers to the subject's phenotype,
#'         obtained after stepwise correlation calculations, are referred to for more details:
#'
#'
#'
#' @examples out<-GE_CV(pheno=pheno,geno=geno,env=env_info,para=envMeanPara,Para_Name="PTT",
#'                        depend="maker",fold=2,reshuffle=5,methods="RM.G",
#'                        ms1=ms1,ms2=ms2)
GP_CV<-function(pheno,geno,env,para,Para_Name,depend=NULL,fold=NULL,reshuffle=NULL,
                  model,methods=NULL,ms1=NULL,ms2=NULL,ENalpha=NULL,GBM_params=NULL,
                nIter=NULL,burnIn=NULL,thin=NULL,SVM_cost=NULL,gamma=NULL,GBM_rounds=NULL){
  if(is.null(depend)){
    depend="norm"
  }
  if(is.null(ENalpha)){
    ENalpha=0.5
  }
  if(is.null(methods)){
    methods="RM.GE"
  }
  if(is.null(fold)){
    fold= 2
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
  if(is.null(nIter)){nIter=5000}
  if(is.null(burnIn)){burnIn=1000}
  if(is.null(thin)){thin=3}
  if(is.null(SVM_cost)){SVM_cost=10}
  if(is.null(gamma)){gamma=0.001}
  if(is.null(GBM_rounds)){GBM_rounds=100L}

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
      return(list(outforfigure,r_within,r_across));
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
      genotype.impute=rrBLUP::A.mat(genotype_1[,-1],max.missing=0.5,impute.method="mean",return.imputed=T);
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
          if(model == "rrBLUP"){
            y0=intercept;
            ans<-rrBLUP::mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
            e=as.matrix(ans$u)
            G.pred=Marker[id.V,]
            y_pred=as.matrix(G.pred) %*% e
            GEBV.inter=c(y_pred[,1])+c(ans$beta);
            y0=slope;
            ans<-rrBLUP::mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
            e=as.matrix(ans$u)
            G.pred=Marker[id.V,]
            y_pred=as.matrix(G.pred) %*% e
            GEBV.slope=c(y_pred[,1])+c(ans$beta);
          }else if(model == "LASSO"){
            y0=intercept;
            cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                      family="gaussian",alpha=1)
            lambda_min <- cv.fit$lambda.min
            coef<-coef(cv.fit)
            e=as.matrix(coef@x[-1])
            G.pred=Marker[id.V,]
            selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            G.pred_LASSO <- G.pred[, selected_features]
            y_pred=as.matrix(G.pred_LASSO) %*% e
            GEBV.inter=c(y_pred[,1])+c(coef@x[1]);

            y0=slope;
            cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                      family="gaussian",alpha=1)
            lambda_min <- cv.fit$lambda.min
            coef<-coef(cv.fit)
            e=as.matrix(coef@x[-1])
            G.pred=Marker[id.V,]
            selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            G.pred_LASSO <- G.pred[, selected_features]
            y_pred=as.matrix(G.pred_LASSO) %*% e
            GEBV.slope=c(y_pred[,1])+c(coef@x[1]);
          }else if(model == "RR"){
            y0=intercept;
            cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                      alpha=0)
            lambda_min <- cv.fit$lambda.min
            coef<-coef(cv.fit)
            e=as.matrix(coef@x[-1])
            G.pred=Marker[id.V,]
            selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            G.pred_RR <- G.pred[, selected_features]
            y_pred=as.matrix(G.pred_RR) %*% e
            GEBV.inter=c(y_pred[,1])+c(coef@x[1]);

            y0=slope;
            cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                      alpha=0)
            lambda_min <- cv.fit$lambda.min
            coef<-coef(cv.fit)
            e=as.matrix(coef@x[-1])
            G.pred=Marker[id.V,]
            selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            G.pred_RR <- G.pred[, selected_features]
            y_pred=as.matrix(G.pred_RR) %*% e
            GEBV.slope=c(y_pred[,1])+c(coef@x[1]);
          }
          else if(model == "EN"){
            y0=intercept;
            cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                      family="gaussian",alpha=ENalpha)
            lambda_min <- cv.fit$lambda.min
            coef<-coef(cv.fit)
            e=as.matrix(coef@x[-1])
            G.pred=Marker[id.V,]
            selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            G.pred_EN <- G.pred[, selected_features]
            y_pred=as.matrix(G.pred_EN) %*% e
            GEBV.inter=c(y_pred[,1])+c(coef@x[1]);

            y0=slope;
            cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                      family="gaussian",alpha=ENalpha)
            lambda_min <- cv.fit$lambda.min
            coef<-coef(cv.fit)
            e=as.matrix(coef@x[-1])
            G.pred=Marker[id.V,]
            selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            G.pred_EN <- G.pred[, selected_features]
            y_pred=as.matrix(G.pred_EN) %*% e
            GEBV.slope=c(y_pred[,1])+c(coef@x[1]);
          } else if(model == "RF"){
            y0=intercept;
            GEBV.inter<-randomForest::randomForest(genotype[id.T,-1], y=y0[id.T],
                                                   xtest=Marker[id.V,])$test$predicted

            y0=slope;
            GEBV.slope<-randomForest::randomForest(genotype[id.T,-1], y=y0[id.T],
                                                   xtest=Marker[id.V,])$test$predicted
          }
          else if(model == "SVM"){
            y0=intercept;
            model.inter<-e1071::svm(genotype[id.T,-1], y=y0[id.T], method="nu-regression",
                                    kernel="radial",cost=SVM_cost,gamma=gamma)
            GEBV.inter<-predict(model.inter,Marker[id.V,])

            y0=slope;
            model.slope<-e1071::svm(genotype[id.T,-1], y=y0[id.T]*10000, method="nu-regression",
                                    kernel="radial",cost=SVM_cost,gamma=gamma,max_iter=10)
            GEBV.slope<-predict(model.slope,Marker[id.V,])/10000
          }
          else if(model == "BA" | model == "BC" | model == "BL" | model == "BRR"
                  | model == "BB" ){
            y0=intercept;
            GENO = genotype[,-1]
            y0[id.V]<-NA
            yNa<-y0
            saveAt = stringi::stri_rand_strings(1, 32, '[a-zA-Z]');
            S0=NULL;weights=NULL;R2=0.5;
            if(model =="BA"){
              nIter=nIter;burnIn=burnIn;thin=thin;
              ETA<-list(list(X=GENO,model='BayesA'))}
            else if (model =="BB"){
              nIter=nIter;burnIn=burnIn;thin=thin;
              ETA<-list(list(X=GENO,model='BayesB'))}
            else if (model =="BC"){
              nIter=nIter;burnIn=burnIn;thin=thin;
              ETA<-list(list(X=GENO,model='BayesC'))}
            else if (model =="BL"){
              nIter=nIter;burnIn=burnIn;thin=thin;
              ETA<-list(list(X=GENO,model='BL'))}
            else if (model =="BRR"){
              nIter=nIter;burnIn=burnIn;thin=thin;
              ETA<-list(list(X=GENO,model='BRR'))}
            model.inter=BGLR::BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,
                       saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2,verbose=F)
            GEBV.inter<-model.inter$yHat[id.V]
            removeSaveAt(saveAt)

            y0=slope;
            GENO = genotype[,-1]
            y0[id.V]<-NA
            yNa<-y0
            saveAt = stringi::stri_rand_strings(1, 32, '[a-zA-Z]');
            S0=NULL;weights=NULL;R2=0.5;
            if(model =="BA"){
              nIter=nIter;burnIn=burnIn;thin=thin;
              ETA<-list(list(X=GENO,model='BayesA'))}
            else if (model =="BB"){
              nIter=nIter;burnIn=burnIn;thin=thin;
              ETA<-list(list(X=GENO,model='BayesB'))}
            else if (model =="BC"){
              nIter=nIter;burnIn=burnIn;thin=thin;
              ETA<-list(list(X=GENO,model='BayesC'))}
            else if (model =="BL"){
              nIter=nIter;burnIn=burnIn;thin=thin;
              ETA<-list(list(X=GENO,model='BL'))}
            else if (model =="BRR"){
              nIter=nIter;burnIn=burnIn;thin=thin;
              ETA<-list(list(X=GENO,model='BRR'))}
            model.slope=BGLR::BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,
                       saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2,verbose=F)
            GEBV.slope<-model.slope$yHat[id.V]
            removeSaveAt(saveAt)
          }
          else if(model =="GBLUP"){
            y0=intercept;
            A<-rrBLUP::A.mat(geno[id.V,-1])
            #return(geno[id.V,-1])
            dataF=data.frame(genoID=rownames(geno[id.V,-1]),intcp=y0[id.V])
            MODEL=rrBLUP::kin.blup(data=dataF,geno="genoID",pheno="intcp", GAUSS=FALSE, K=A,
                                   PEV=TRUE,n.core=1,theta.seq=NULL)
            GEBV.inter<-MODEL$pred


            y0=slope;
            A<-rrBLUP::A.mat(geno[id.V,-1])
            dataF=data.frame(genoID=rownames(geno[id.V,-1]),slope=y0[id.V])
            MODEL=rrBLUP::kin.blup(data=dataF,geno="genoID",pheno="slope", GAUSS=FALSE, K=A,
                                   PEV=TRUE,n.core=1,theta.seq=NULL)
            GEBV.slope<-MODEL$pred

          }
          else if(model =="RKHS" | model =="MKRKHS"){
            y0=intercept;
            GENO = genotype[,-1]
            y0[id.V]<-NA
            yNa<-y0

            X = GENO
            p = ncol(X)

            X <- scale(X,center=TRUE,scale=TRUE)
            D <- (as.matrix(dist(X,method='euclidean'))^2)/p
            h <- 0.5
            K=exp(-h*D)
            saveAt=stringi::stri_rand_strings(1, 32, '[a-zA-Z]');
            if(model =="RKHS"){
            ETA= list(list(K=K, model='RKHS'))
            MODEL <- BGLR::BGLR(y=yNa,ETA=ETA,nIter=nIter, burnIn=burnIn, saveAt = saveAt,verbose=F)}
            else if(model =="MKRKHS"){
              h <- sqrt(c(0.2,0.5,0.8))
              #3# Kernel Averaging using BGLR
              ETA <- list(list(K=exp(-h[1]*D),model='RKHS'),
                          list(K=exp(-h[2]*D),model='RKHS'),
                          list(K=exp(-h[3]*D),model='RKHS'))
              MODEL <- BGLR::BGLR(y=yNa,ETA=ETA,nIter=nIter, burnIn=burnIn, saveAt = saveAt,verbose=F)
            }
            GEBV.inter<-MODEL$yHat[id.V]

            y0=slope;
            GENO = genotype[,-1]
            y0[id.V]<-NA
            yNa<-y0

            X = GENO
            p = ncol(X)

            X <- scale(X,center=TRUE,scale=TRUE)
            D <- (as.matrix(dist(X,method='euclidean'))^2)/p
            h <- 0.5
            K=exp(-h*D)
            saveAt=stringi::stri_rand_strings(1, 32, '[a-zA-Z]');

            if(model =="RKHS"){
            ETA= list(list(K=K, model='RKHS'))
            MODEL <- BGLR::BGLR(y=yNa,ETA=ETA,nIter=nIter, burnIn=burnIn, saveAt = saveAt,verbose=F)
          }
          else if(model =="MKRKHS"){
            h <- sqrt(c(0.2,0.5,0.8))
            #3# Kernel Averaging using BGLR
            ETA <- list(list(K=exp(-h[1]*D),model='RKHS'),
                        list(K=exp(-h[2]*D),model='RKHS'),
                        list(K=exp(-h[3]*D),model='RKHS'))
            MODEL <- BGLR::BGLR(y=yNa,ETA=ETA,nIter=nIter, burnIn=burnIn, saveAt = saveAt,verbose=F)
          }
            GEBV.slope<-MODEL$yHat[id.V]
            #return(list(GEBV.inter,GEBV.slope))
            removeSaveAt(saveAt)
          }
          if(model == "GBM"){
            if(is.null(GBM_params)){
            params <- list(boosting="gbdt",objective = "regression",metric = "RMSE",min_data = 1L,
                           learning_rate = 0.01,num_iterations=1000,num_leaves=3,max_depth=-1,
                           early_stopping_round=50L,cat_l2=10,skip_drop=0.5,drop_rate=0.5,
                           cat_smooth=5)
            }
            #print(params)
            y0=intercept;
            dtrain <- lightgbm::lgb.Dataset(genotype[id.T,-1], label = y0[id.T])
            dtest <- lightgbm::lgb.Dataset.create.valid(dtrain, genotype[id.V,-1], label = y0[id.V])
            valids <- list(test = dtest)
            model.inter <- lightgbm::lgb.train(params = params,data = dtrain
                                               ,valids = valids,nrounds = GBM_rounds,verbose = -1)
            GEBV.inter <- predict(model.inter, Marker[id.V,])

            y0=slope;
            dtrain <- lightgbm::lgb.Dataset(genotype[id.T,-1], label = y0[id.T])
            dtest <- lightgbm::lgb.Dataset.create.valid(dtrain, genotype[id.V,-1], label = y0[id.V])
            valids <- list(test = dtest)
            model.slope <- lightgbm::lgb.train(params = params,data = dtrain,
                                               valids = valids,nrounds = GBM_rounds,verbose = -1)
            GEBV.slope <- predict(model.slope, Marker[id.V,])
            #rm(ls(model.slope,model.slope))
            #return(model.inter)
          }
          yhat.para=matrix(999,length(id.V),n.para);
          yobs.para=matrix(999,length(id.V),n.para);

          for(j in 1:n.para)
          {
            yhat=GEBV.inter+GEBV.slope*para[j,enp];
            yobs=pheno[id.V,j+1];
            #return(list(yhat,yhat.para[,1]))
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
      genotype.impute=rrBLUP::A.mat(genotype_1[,-1],max.missing=0.5,impute.method="mean",return.imputed=T);
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

            if(model == "rrBLUP"){
            y0=intercept;
            ans<-rrBLUP::mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
            e=as.matrix(ans$u)
            G.pred=Marker[id.V,]
            y_pred=as.matrix(G.pred) %*% e
            GEBV.inter=c(y_pred[,1])+c(ans$beta);

            y0=slope;
            ans<-rrBLUP::mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
            e=as.matrix(ans$u)
            G.pred=Marker[id.V,]
            y_pred=as.matrix(G.pred) %*% e
            GEBV.slope=c(y_pred[,1])+c(ans$beta);

            }
            else if (model == "LASSO"){
              y0=intercept;
              cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                        family="gaussian",alpha=1)
              lambda_min <- cv.fit$lambda.min
              coef<-coef(cv.fit)
              e=as.matrix(coef@x[-1])
              G.pred=Marker[id.V,]
              selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
              G.pred_LASSO <- G.pred[, selected_features]
              y_pred=as.matrix(G.pred_LASSO) %*% e
              GEBV.inter=c(y_pred[,1])+c(coef@x[1]);

              y0=slope;
              cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                        family="gaussian",alpha=1)
              lambda_min <- cv.fit$lambda.min
              coef<-coef(cv.fit)
              e=as.matrix(coef@x[-1])
              G.pred=Marker[id.V,]
              selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
              G.pred_LASSO <- G.pred[, selected_features]
              y_pred=as.matrix(G.pred_LASSO) %*% e
              GEBV.slope=c(y_pred[,1])+c(coef@x[1]);
            }
            else if (model == "RR"){
              y0=intercept;
              cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                        alpha=0)
              lambda_min <- cv.fit$lambda.min
              coef<-coef(cv.fit)
              e=as.matrix(coef@x[-1])
              G.pred=Marker[id.V,]
              selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
              G.pred_LASSO <- G.pred[, selected_features]
              y_pred=as.matrix(G.pred_LASSO) %*% e
              GEBV.inter=c(y_pred[,1])+c(coef@x[1]);

              y0=slope;
              cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                        alpha=0)
              lambda_min <- cv.fit$lambda.min
              coef<-coef(cv.fit)
              e=as.matrix(coef@x[-1])
              G.pred=Marker[id.V,]
              selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
              G.pred_LASSO <- G.pred[, selected_features]
              y_pred=as.matrix(G.pred_LASSO) %*% e
              GEBV.slope=c(y_pred[,1])+c(coef@x[1]);
            }
            else if (model == "EN"){
              y0=intercept;
              cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                        family="gaussian",alpha=ENalpha)
              lambda_min <- cv.fit$lambda.min
              coef<-coef(cv.fit)
              e=as.matrix(coef@x[-1])
              G.pred=Marker[id.V,]
              selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
              G.pred_EN <- G.pred[, selected_features]
              y_pred=as.matrix(G.pred_EN) %*% e
              GEBV.inter=c(y_pred[,1])+c(coef@x[1]);

              y0=slope;
              cv.fit<-glmnet::cv.glmnet(y=y0[id.T],x=genotype[id.T,-1],
                                        family="gaussian",alpha=ENalpha)
              lambda_min <- cv.fit$lambda.min
              coef<-coef(cv.fit)
              e=as.matrix(coef@x[-1])
              G.pred=Marker[id.V,]
              selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
              G.pred_EN <- G.pred[, selected_features]
              y_pred=as.matrix(G.pred_EN) %*% e
              GEBV.slope=c(y_pred[,1])+c(coef@x[1]);
            }
            else if (model == "SVM"){
              y0=intercept;
              model.inter<-e1071::svm(genotype[id.T,-1], y=y0[id.T], method="nu-regression",
                                      kernel="radial",cost=10,gamma=0.001)
              GEBV.inter<-predict(model.inter,Marker[id.V,])

              y0=slope;
              model.slope<-e1071::svm(genotype[id.T,-1], y=y0[id.T], method="nu-regression",
                                      kernel="radial",cost=10,gamma=0.001)
              GEBV.slope<-predict(model.slope,Marker[id.V,])
            }
            else if (model == "RF"){
              y0=intercept;
              GEBV.inter<-randomForest::randomForest(genotype[id.T,-1], y=y0[id.T],
                                                     xtest=Marker[id.V,])$test$predicted

              y0=slope;
              GEBV.slope<-randomForest::randomForest(genotype[id.T,-1], y=y0[id.T],
                                                     xtest=Marker[id.V,])$test$predicted
            }
            else if(model =="GBLUP"){
              y0=intercept;
              A<-rrBLUP::A.mat(geno[id.V,-1])
              dataF=data.frame(genoID=rownames(geno[id.V,-1]),intcp=y0[id.V])
              MODEL=rrBLUP::kin.blup(data=dataF,geno="genoID",pheno="intcp", GAUSS=FALSE, K=A,
                                     PEV=TRUE,n.core=1,theta.seq=NULL)
              GEBV.inter<-MODEL$pred


              y0=slope;
              A<-rrBLUP::A.mat(geno[id.V,-1])

              dataF=data.frame(genoID=rownames(geno[id.V,-1]),slope=y0[id.V])
              MODEL=rrBLUP::kin.blup(data=dataF,geno="genoID",pheno="slope", GAUSS=FALSE, K=A,
                                     PEV=TRUE,n.core=1,theta.seq=NULL)
              GEBV.slope<-MODEL$pred
            }
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
  }else if(depend =="Marker"){
    genotype.match=match(pheno[,1],geno[,1])
    genotype_1=geno[genotype.match,];
    genotype.impute=rrBLUP::A.mat(genotype_1[,-1],max.missing=0.5,
                                  impute.method="mean",return.imputed=T);
    SFI=cbind(genotype_1[,1],genotype.impute$imputed);
    genotype=matrix(suppressWarnings(as.numeric(SFI)),nrow=nrow(SFI))
    Marker=genotype[,-1];
    n.marker=dim(Marker)[2];
    if(methods=="RM.E"){
      effect=matrix(999,n.marker,n.para)
      intercept=matrix(999,1,n.para)
      for(i in 1:n.para)
      {
        fit=rrBLUP::mixed.solve(pheno[,1+i],Z=Marker)
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
      cor.within<-matrix(999,reshuffle,n.para)
      cor.all<-numeric()
      if(model == "rrBLUP" | model == "RR"){
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
            if(model == "rrBLUP"){
            fit=rrBLUP::mixed.solve(pheno[id.T,1+k],Z=Marker[id.T,])
            effect[,k]=fit$u
            intercept[,k]=fit$beta
            }
            if(model == "RR"){
              cv.fit<-glmnet::cv.glmnet(y=pheno[id.T,1+k],x=Marker[id.T,],
                                        alpha=0)
              lambda_min <- cv.fit$lambda.min
              coef<-coef(cv.fit)
              effect[,k]=coef@x[-1]
              intercept[,k]=coef@x[1]
            }
          }
          ##Slope###
          effect.hat=numeric();
          for(j in 1:n.marker)
          {
            x1=para[,enp]
            y1=effect[j,]
            coe=lm(y~x,data=data.frame(x=x1,y=y1))
            y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*para[,enp]
            effect.hat=rbind(effect.hat,y.hat)
          }

          ##Environment intercept####
          reg.intercept=data.frame(x=as.numeric(para[,enp]),y=as.vector(intercept));
          coe.intercept=lm(y~x,data=reg.intercept)
          y.intercept=summary(coe.intercept)$coefficients[1,1]+
            summary(coe.intercept)$coefficients[2,1]*as.numeric(para[,enp]);

          ###All the predicted slope and intercept
          yhat.para=matrix(999,length(id.V),n.para)
          yobs.para=matrix(999,length(id.V),n.para)
          for(j in 1:n.para)
          {
            yobs=pheno[id.V,j+1]
            yhat=y.intercept[j]+(as.matrix(Marker[id.V,])%*%as.matrix(effect.hat))[,j]
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
        r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1]
        #Correlation across environment 50 times
        r_across=mean(cor.all)
        #Observation and prediction last time
        colnames(cor.within)=colnames(pheno)[-1]
        r_within=cor.within
        r_across=cor.all
        #return(list(r_within,r_across))
        #return(list(outforfigure,r_within,r_across))
        outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                                pre=as.vector(yhat.whole.cross),
                                col=rep(coloo,times=rep(n.line,n.para)),
                                para=rep(colnames(pheno)[-1],times=rep(n.line,n.para)))
        return(list(outforfigure,r_within,r_across))
      }
      else{
        cor.within<-matrix(999,reshuffle,n.para)
        cor.all<-numeric()
        for(r in 1:reshuffle){
          cross = sample(rep(1:fold,each=ceiling(n.line/fold)),n.line)
          yhat.whole.cross<-data.frame()
          yobs.whole.cross<-data.frame()
          for(f in 1:fold)
          {
            id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f]
            ##marker effect###
            #Test SNP matrix
            Test<-as.matrix(geno[id.T,-1])
            n <- nrow(Test)
            m <- ncol(Test)
            new_n <- n * n.para
            new_m <- m * n.para
            M <- matrix(0, nrow = new_n, ncol = new_m)
            for (i in 1:n.para) {
              row_start <- (i - 1) * n + 1
              row_end <- i * n
              col_start <- (i - 1) * m + 1
              col_end <- i * m
              M[row_start:row_end, col_start:col_end] <- Test
            }
            #Pred SNP matrix
            Pred<-as.matrix(geno[id.V,-1])
            n <- nrow(Pred)
            m <- ncol(Pred)
            new_n <- n * n.para
            new_m <- m * n.para
            N <- matrix(0, nrow = new_n, ncol = new_m)
            for (i in 1:n.para) {
              row_start <- (i - 1) * n + 1
              row_end <- i * n
              col_start <- (i - 1) * m + 1
              col_end <- i * m
              N[row_start:row_end, col_start:col_end] <- Pred
            }

            n <- nrow(as.matrix(geno[,-1]))
            m <- ncol(as.matrix(geno[,-1]))
            GENO <- matrix(0, nrow = n * n.para, ncol = m * n.para)
            for (i in 1:n.para) {
              row_start <- (i - 1) * n + 1
              row_end <- i * n
              col_start <- (i - 1) * m + 1
              col_end <- i * m
              GENO[row_start:row_end, col_start:col_end] <- as.matrix(geno[,-1])
            }

            #pheno1
            yNa<-pheno[which(as.character(pheno$line_code)%in%
                                  c("line_code",as.character(geno[id.T,1]))),]
            yNas <- yNa[order(factor(yNa$line_code, levels = geno[id.T,1])), ]
            yNal <- yNas %>%tidyr::pivot_longer(cols = -line_code,
                                                         names_to = "Envs",
                                                         values_to = "Obs")
            yNal<-as.data.frame(yNal)
            yNa<-yNal[order(yNal[,2]), ]
            #pheno2
            y_Prd <- pheno[which(as.character(pheno$line_code)%in%
                                   c("line_code",as.character(geno[id.V,1]))),]
            y_Prds <- y_Prd[order(factor(y_Prd$line_code, levels = geno[id.V,1])), ]
            y_Prdl <- y_Prds %>%tidyr::pivot_longer(cols = -line_code,
                                                         names_to = "Envs",
                                                         values_to = "Obs")
            y_Prdl<-as.data.frame(y_Prdl)
            y_Prd<-y_Prdl[order(y_Prdl[,2]), ]
            #pheno3
            y <- pheno[which(as.character(pheno$line_code)%in%
                                   c("line_code",as.character(geno[,1]))),]
            ys <- y[order(factor(y$line_code, levels = geno[,1])), ]
            yl <- ys %>%tidyr::pivot_longer(cols = -line_code,
                                                    names_to = "Envs",
                                                    values_to = "Obs")
            yl<-as.data.frame(yl)
            y<-yl[order(yl[,2]),]

            if (model == "RF"){
              Prd<-randomForest::randomForest(M, y=yNa[,ncol(yNa)],xtest=N)$test$predicted
              PrdM<-data.frame(y_Prd, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = line_code,names_from = Envs,values_from = Prd)
              PrdM_wide<-as.data.frame(PrdM_wide)
            }
            else if (model == "LASSO"){
              cv.fit <- glmnet::cv.glmnet(M,yNa[,ncol(yNa)],family="gaussian",alpha=1)
              lambda_min <- cv.fit$lambda.min
              coef<-coef(cv.fit)
              e=as.matrix(coef@x[-1])
              G.pred=N
              selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)
              G.pred_LASSO <- G.pred[, selected_features]
              y_pred=as.matrix(G.pred_LASSO) %*% e
              Prd=c(y_pred[,1])+c(coef@x[1])
              PrdM<-data.frame(y_Prd, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(
                id_cols = line_code,
                names_from = Envs,
                values_from = Prd
              )
              PrdM_wide<-as.data.frame(PrdM_wide)
            }
            else if (model == "EN"){
              cv.fit <- glmnet::cv.glmnet(M,yNa[,ncol(yNa)],family="gaussian",alpha=ENalpha)
              lambda_min <- cv.fit$lambda.min
              coef<-coef(cv.fit)
              e=as.matrix(coef@x[-1])
              G.pred=N
              selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
              G.pred_EN <- G.pred[, selected_features]
              y_pred=as.matrix(G.pred_LASSO) %*% e
              Prd=c(y_pred[,1])+c(coef@x[1])
              PrdM<-data.frame(y_Prd, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(
                id_cols = line_code,
                names_from = Envs,
                values_from = Prd
              )
              PrdM_wide<-as.data.frame(PrdM_wide)
            }
            else if (model == "GBLUP"){
              GENO<-as.data.frame(GENO)
              rowname<-paste(y[,1],y[,2], sep="-")
              rownames(GENO)<-rowname
              #A<-rrBLUP::A.mat(GENO)
              y2<-y %>%
                group_by(Envs) %>%
                mutate(row_number = row_number()) %>%
                mutate(Obs = ifelse(row_number %in% id.V, NA, Obs)) %>%
                select(-row_number)
              y2<-as.data.frame(y2)
              dataF=data.frame(genoID=rownames(GENO),yield=y2[,ncol(y2)])
              MODEL<-rrBLUP::kin.blup(data=dataF,geno="genoID",pheno="yield", GAUSS=FALSE, K=A,
                                     PEV=TRUE,n.core=1,theta.seq=NULL)
              Prd<-MODEL$pred
              PrdM<-data.frame(y2, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(
                id_cols = line_code,
                names_from = Envs,
                values_from = Prd
              )
              PrdM_wide<-as.data.frame(PrdM_wide)
              PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% geno[id.V,1], ]
              PrdM_wide<-as.data.frame(PrdM_wide)
            }
            else if (model == "BA" | model == "BC" | model == "BL" | model == "BRR"
               | model == "BB" ){
              saveAt = stringi::stri_rand_strings(1, 32, '[a-zA-Z]');
              S0=NULL;weights=NULL;R2=0.5;
              if(model == "BA"){
                nIter=nIter;burnIn=burnIn;thin=thin;
                ETA<-list(list(X=GENO,model='BayesA'))}
              else if (model == "BB"){
                nIter=nIter;burnIn=burnIn;thin=thin;
                ETA<-list(list(X=GENO,model='BayesB'))}
              else if (model == "BC"){
                nIter=nIter;burnIn=burnIn;thin=thin;
                ETA<-list(list(X=GENO,model='BayesC'))}
              else if (model == "BL"){
                nIter=nIter;burnIn=burnIn;thin=thin;
                ETA<-list(list(X=GENO,model='BL'))}
              else if (model == "BRR"){
                nIter=nIter;burnIn=burnIn;thin=thin;
                ETA<-list(list(X=GENO,model='BRR'))}
              y2<-y %>%
                group_by(Envs) %>%
                mutate(row_number = row_number()) %>%
                mutate(Obs = ifelse(row_number %in% id.V, NA, Obs)) %>%
                select(-row_number)
              y2<-as.data.frame(y2)
              model.inter=BGLR::BGLR(y=y2[,ncol(y2)],ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,
                                     saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2,verbose=F)
              #return(list(y,model.inter))
              PrdM<-data.frame(y2, Prd = as.numeric(model.inter$yHat))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(
                id_cols = line_code,
                names_from = Envs,
                values_from = Prd
              )
              PrdM_wide<-as.data.frame(PrdM_wide)
              PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% geno[id.V,1], ]
              PrdM_wide<-as.data.frame(PrdM_wide)
              removeSaveAt(saveAt)
            }
            else if (model == "RKHS" | model == "MKRKHS"){
              GENO<-as.data.frame(GENO)
              rowname<-paste(y[,1],y[,2], sep="-")
              rownames(GENO)<-rowname
              y2<-y %>%
                group_by(Envs) %>%
                mutate(row_number = row_number()) %>%
                mutate(Obs = ifelse(row_number %in% id.V, NA, Obs)) %>%
                select(-row_number)
              y2<-as.data.frame(y2)
              X <- GENO
              p <- ncol(X)
              X <- scale(X,center=TRUE,scale=TRUE)
              D <- (as.matrix(dist(X,method='euclidean'))^2)/p
              h <- 0.5
              K <- exp(-h*D)
              saveAt=stringi::stri_rand_strings(1, 32, '[a-zA-Z]')
              if(model =="RKHS"){
                ETA= list(list(K=K, model='RKHS'))
                model_RKHS<-BGLR::BGLR(y=y2[,ncol(y2)],ETA=ETA,nIter=nIter, burnIn=burnIn,
                                       saveAt = saveAt,verbose=F)
                }
              else if(model =="MKRKHS"){
                h <- sqrt(c(0.2,0.5,0.8))
                #3# Kernel Averaging using BGLR
                ETA <- list(list(K=exp(-h[1]*D),model='RKHS'),
                            list(K=exp(-h[2]*D),model='RKHS'),
                            list(K=exp(-h[3]*D),model='RKHS'))
                model_RKHS<-BGLR::BGLR(y=y2[,ncol(y2)],ETA=ETA,nIter=nIter, burnIn=burnIn,
                                       saveAt = saveAt,verbose=F)
              }
              PrdM<-data.frame(y, Prd = as.numeric(model_RKHS$yHat))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(
                id_cols = line_code,
                names_from = Envs,
                values_from = Prd
              )
              PrdM_wide<-as.data.frame(PrdM_wide)
              PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% geno[id.V,1], ]
              PrdM_wide<-as.data.frame(PrdM_wide)
              removeSaveAt(saveAt)
            }
            else if (model == "SVM"){
              model_SVM <- e1071::svm(M,yNa[,ncol(yNa)],method="nu-regression",
                                  kernel="radial",cost=10,gamma=0.001)
              Prd <- predict(model_SVM,N)
              PrdM<-data.frame(y_Prd, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = line_code,names_from = Envs,values_from = Prd)
              PrdM_wide<-as.data.frame(PrdM_wide)
            }
            else if (model == "GBM"){
              if(is.null(GBM_params)){
                params <- list(boosting="gbdt",objective = "regression",metric = "RMSE",min_data = 1L,
                               learning_rate = 0.01,num_iterations=1000,num_leaves=3,max_depth=-1,
                               early_stopping_round=50L,cat_l2=10,skip_drop=0.5,drop_rate=0.5,
                               cat_smooth=5)
              }
              dtrain <- lightgbm::lgb.Dataset(M, label = yNa[,ncol(yNa)])
              dtest <- lightgbm::lgb.Dataset.create.valid(dtrain, N, label = y_Prd[,ncol(y_Prd)])
              valids <- list(test = dtest)
              model_GBM <- lightgbm::lgb.train(params = params,data = dtrain,
                                           valids = valids,nrounds = GBM_rounds,verbose = -1)
              Prd <- predict(model_GBM, N)
              PrdM<-data.frame(y_Prd, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = line_code,names_from = Envs,values_from = Prd)
              PrdM_wide<-as.data.frame(PrdM_wide)
            }
            yhat.whole.cross<-rbind(yhat.whole.cross,PrdM_wide)
            #print(dim(yhat.whole.cross))
            #return(yhat.whole.cross)
          }
          #need improve
          yhat.whole.cross <- rbind(yhat.whole.cross[order(factor(yhat.whole.cross$line_code, levels = ys[,1])), ],yNas)
          yhat.whole.cross <- aggregate(. ~ line_code, data = yhat.whole.cross, FUN = mean)
          yobs.whole.cross <- ys
          yhat.whole.cross <- yhat.whole.cross[order(factor(yhat.whole.cross$line_code, levels = yobs.whole.cross[,1])), ]
          #return(list(yobs.whole.cross,yhat.whole.cross))

          yhat.whole.crossl<-tidyr::gather(yhat.whole.cross, key = "Envs", value = "Prd", -line_code)
          yobs.whole.crossl<-tidyr::gather(yobs.whole.cross, key = "Envs", value = "Obs", -line_code)
          cor.all<-c(cor.all,cor(as.vector(yhat.whole.crossl[,3]),
                                as.vector(yobs.whole.crossl[,3]),
                                use = "complete.obs"))
          #print(cor.all)
          for(k in 1:n.para){
            cor.within[r,k]<-cor(yhat.whole.cross[,k+1],yobs.whole.cross[,k+1],use = "complete.obs")
          }
          #return(cor.within)
        }
      #Correlation within environment 50 times
      r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1]
      #Correlation across environment 50 times
      r_across=mean(cor.all)
      #Observation and prediction last time
      colnames(cor.within)=colnames(pheno)[-1]
      r_within=cor.within
      r_across=cor.all
      outforfigure<-rbind(yobs.whole.crossl,yobs.whole.crossl)
      outforfigure<-data.frame(outforfigure,col=rep(coloo,times=rep(n.line,n.para)))
      return(list(outforfigure,r_within,r_across))
      #return(list(outforfigure,r_within,r_across))
      outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                              pre=as.vector(yhat.whole.cross),
                              col=rep(coloo,times=rep(n.line,n.para)),
                              para=rep(colnames(pheno)[-1],times=rep(n.line,n.para)))
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
            if(model == "rrBLUP"){
              fit=rrBLUP::mixed.solve(pheno[id.T,1+k],Z=Marker[id.T,])
              effect[,k]=fit$u
              intercept[,k]=fit$beta}
            if(model == "RR"){
              cv.fit<-glmnet::cv.glmnet(y=pheno[id.T,1+k],x=Marker[id.T,],
                                        alpha=0)
              lambda_min <- cv.fit$lambda.min
              coef<-coef(cv.fit)
              effect[,k]=coef@x[-1]
              intercept[,k]=coef@x[1]
            }
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
            yhat=y.intercept+(as.matrix(Marker[id.V,])%*%as.matrix(effect.hat))[,j];
            yobs=pheno[id.V,kk+1];
            obs.para=cbind(obs.para,yobs);
            pre.para=cbind(pre.para,yhat);
          }
          #########

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
}
removeSaveAt <- function(saveAt) {
  Sys.sleep(2)
  files <- list.files('.', paste0("^", saveAt, '.*(lambda|mu|varE|varB|varU|ScaleBayesA|parBayesB|parBayesC)\\.dat$'))
  tryCatch({
    invisible(lapply(files, file.remove))
  }, error = function(e) {

  })
}
