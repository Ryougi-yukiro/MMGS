#' Cross_forecast
#'
<<<<<<< HEAD
#' @param pheno_train  pheno data frame used for training
#' @param depend use which depend
#' @param predict.method use which predict model, you can see this in MMGP
#' @param geno geno martrix
#' @param env_info env info
#' @param predict_env env need predict
#' @param reshuffle times for duplication
=======
#' @param geno_train
#' @param pheno_train
#' @param envs
#' @param depend
#' @param predict.method
>>>>>>> 23b0efae316a77181b94cb5eff39d6257393dcf5
#'
#' @return predict value
#' @export
#'
#' @examples
<<<<<<< HEAD
#' \dontrun{
#' Cross_forecast(geno,pheno_train,env_info,predict_env,reshuffle,depend,predict.method)
#' }
Cross_forecast<-function(geno,pheno_train,env_info,predict_env,reshuffle,depend,predict.method){
=======
Cross_forecast<-function(geno,pheno_train,env_info,predict_env,reshuffle,depend="Marker",predict.method = "GBLUP"){
>>>>>>> 23b0efae316a77181b94cb5eff39d6257393dcf5
  if(is.null(nIter)){nIter=5000}
  if(is.null(burnIn)){burnIn=1000}
  if(is.null(thin)){thin=3}
  if(is.null(SVM_cost)){SVM_cost=10}
  if(is.null(gamma)){gamma=0.001}
  if(is.null(GBM_rounds)){GBM_rounds=100L}
  if(is.null(ENalpha)){ENalpha=0.5}
  if(is.null(reshuffle)){reshuffle= 5}
  geno<-geno_train
  envs<-env_info
  n.para<-dim(envs)[1]
  Out<-data.frame()
for (i in reshuffle){
  if(depend == "Marker"){
    #GENO
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
    #pheno
    y <- pheno[which(as.character(pheno$line_code)%in%
                       c("line_code",as.character(geno[,1]))),]
    ys <- y[order(factor(y$line_code, levels = geno[,1])), ]
    yl <- ys %>%tidyr::pivot_longer(cols = -line_code,
                                    names_to = "Envs",
                                    values_to = "Obs")
    yl<-as.data.frame(yl)
    y<-yl[order(yl[,2]),]
    if (model == "RF"){
      Prd<-randomForest::randomForest(GENO, y=y[,ncol(y)],xtest=GENO)$test$predicted
      PrdM<-data.frame(y, Prd = as.numeric(Prd))
      PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = line_code,names_from = Envs,values_from = Prd)
      PrdM_wide<-as.data.frame(PrdM_wide)
    }
    else if (model == "LASSO"){
      cv.fit <- glmnet::cv.glmnet(GENO, y=y[,ncol(y)],family="gaussian",alpha=1)
      lambda_min <- cv.fit$lambda.min
      coef<-coef(cv.fit)
      e<-as.matrix(coef@x[-1])
      G.pred<-GENO
      selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)
      G.pred_LASSO <- G.pred[, selected_features]
      y_pred=as.matrix(G.pred_LASSO) %*% e
      Prd=c(y_pred[,1])+c(coef@x[1])
      PrdM<-data.frame(y, Prd = as.numeric(Prd))
      PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = line_code,names_from = Envs,values_from = Prd)
      PrdM_wide<-as.data.frame(PrdM_wide)
    }
    else if (model == "EN"){
      cv.fit <- glmnet::cv.glmnet(GENO, y=y[,ncol(y)],family="gaussian",alpha=ENalpha)
      lambda_min <- cv.fit$lambda.min
      coef<-coef(cv.fit)
      e<-as.matrix(coef@x[-1])
      G.pred<-GENO
      selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)
      G.pred_EN <- G.pred[, selected_features]
      y_pred<-as.matrix(G.pred_LASSO) %*% e
      Prd<-c(y_pred[,1])+c(coef@x[1])
      PrdM<-data.frame(y, Prd = as.numeric(Prd))
      PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = line_code,names_from = Envs,values_from = Prd)
      PrdM_wide<-as.data.frame(PrdM_wide)
    }
    else if (model == "GBLUP"){
      GENO<-as.data.frame(GENO)
      rowname<-paste(y[,1],y[,2], sep="-")
      rownames(GENO)<-rowname
      dataF=data.frame(genoID=rownames(GENO),yield=y[,ncol(y)])
      MODEL<-rrBLUP::kin.blup(data=dataF,geno="genoID",pheno="yield", GAUSS=FALSE, K=A,
                              PEV=TRUE,n.core=1,theta.seq=NULL)
      Prd<-MODEL$pred
      PrdM<-data.frame(y, Prd = as.numeric(Prd))
      PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = line_code,names_from = Envs,values_from = Prd)
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
      model.inter=BGLR::BGLR(y=y[,ncol(y)],ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,
                             saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2,verbose=F)
      PrdM<-data.frame(y, Prd = as.numeric(model.inter$yHat))
      PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = line_code,names_from = Envs,values_from = Prd)
      PrdM_wide<-as.data.frame(PrdM_wide)
      removeSaveAt(saveAt)
    }
    else if (model == "RKHS" | model == "MKRKHS"){
      GENO<-as.data.frame(GENO)
      rowname<-paste(y[,1],y[,2], sep="-")
      rownames(GENO)<-rowname
      X <- GENO
      p <- ncol(X)
      X <- scale(X,center=TRUE,scale=TRUE)
      D <- (as.matrix(dist(X,method='euclidean'))^2)/p
      h <- 0.5
      K <- exp(-h*D)
      saveAt=stringi::stri_rand_strings(1, 32, '[a-zA-Z]')
      if(model =="RKHS"){
        ETA= list(list(K=K, model='RKHS'))
        model_RKHS<-BGLR::BGLR(y=y[,ncol(y)],ETA=ETA,nIter=nIter, burnIn=burnIn,
                               saveAt = saveAt,verbose=F)
      }
      else if(model =="MKRKHS"){
        h <- sqrt(c(0.2,0.5,0.8))
        #3# Kernel Averaging using BGLR
        ETA <- list(list(K=exp(-h[1]*D),model='RKHS'),
                    list(K=exp(-h[2]*D),model='RKHS'),
                    list(K=exp(-h[3]*D),model='RKHS'))
        model_RKHS<-BGLR::BGLR(y=y[,ncol(y)],ETA=ETA,nIter=nIter, burnIn=burnIn,
                               saveAt = saveAt,verbose=F)
      }
      PrdM<-data.frame(y, Prd = as.numeric(model_RKHS$yHat))
      PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = line_code,names_from = Envs,values_from = Prd)
      PrdM_wide<-as.data.frame(PrdM_wide)
      removeSaveAt(saveAt)
    }
    else if (model == "SVM"){
      model_SVM <- e1071::svm(GENO,y[,ncol(y)],method="nu-regression",
                              kernel="radial",cost=10,gamma=0.001)
      Prd <- predict(model_SVM,GENO)
      PrdM<-data.frame(y, Prd = as.numeric(Prd))
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
      dtrain <- lightgbm::lgb.Dataset(GENO, label = y[,ncol(y)])
      dtest <- lightgbm::lgb.Dataset.create.valid(dtrain, GENO, label = y[,ncol(y)])
      valids <- list(test = dtest)
      model_GBM <- lightgbm::lgb.train(params = params,data = dtrain,
                                       valids = valids,nrounds = GBM_rounds,verbose = -1)
      Prd <- predict(model_GBM, GENO)
      PrdM<-data.frame(y, Prd = as.numeric(Prd))
      PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = line_code,names_from = Envs,values_from = Prd)
      PrdM_wide<-as.data.frame(PrdM_wide)
    }
  }
  if(depend == "Norm"){

  }
  #统计多次重复取平均结果
  Prd_Env<-PrdM_wide[[predict_env]]
  Out<-rbind(Out,Prd_Env)
<<<<<<< HEAD
}
  Out<-rbind(Out,PrdM_wide)
  Out<-stats::aggregate(. ~ line_code, data = Out, FUN = mean)
=======
  }
>>>>>>> 23b0efae316a77181b94cb5eff39d6257393dcf5
}
