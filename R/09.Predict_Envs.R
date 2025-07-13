#' @title MMPrdM
#'
#' @description
#' MMPrdM (Multi-environment Multi-marker Prediction Model) is a function
#' for genomic prediction that includes cross-environment prediction.
#' It leverages genetic and environmental information to predict phenotype.
#' The function supports various genomic prediction methods and models,
#' providing flexibility in choosing the most suitable approach
#' based on the data characteristics and user preferences.
#'
#' @param geno Matrix (n x m) of genotypes for the training population:
#'             n lines with m markers.
#'             Genotypes should be coded -1, 0, 1.
#'             Missing data are not allowed.
#' @param pheno Vector (n x j) of "phenotypes", i.e. observations or pre-processed, corrected values.
#' @param env Data.frame (j x l) of "environmental information", i.e.
#' @param para Data.frame returned from envMeanPara function.
#' @param Para_Name The most relevant environmental silvers to the subject's phenotype,
#'                  obtained after stepwise correlation calculations, are referred to for more details:
#' @param depend The options for genomic breeding within different environment.The available options are:
#'               1.Reaction Norm:
#'               2.PEI:
#' @param reshuffle Number of independent replicates for the Predict.
#'                  Smallest value recommended is reshuffle = 3.
#' @param methods RM.G RM.GE RM.E
#' @param model The options for genomic breeding cross vadiation methods. The available options are:<br>
#' 1.GBLUP: performs G-BLUP using a marker-based relationship matrix, implemented through BGLR R-library.<br>
#'              Equivalent to ridge regression (RR-BLUP) of marker effects.<br>
#' 2.RR: ridge regression, using package glmnet.
#'              In theory, strictly equivalent to gblup.<br>
#' 3.LASSO: Least Absolute Shrinkage and Selection Operator is another penalized regression methods<br>
#'              which yield more shrinked estimates than RR.Run by glmnet library.<br>
#'              4.EN: Elastic Net (Zou and Hastie, 2005),
#'              which is a weighted combination of RR and LASSO, using glmnet library<br>
#'              5.rrBLUP:the RR-BLUP mixed model(Endelman ,2011). <br>
#'              One application is to estimate marker effects by ridge regression;<br>
#'              alternatively, BLUPs can be calculated based on an <br>
#'              additive relationship matrix or a Gaussian kernel.<br>
#'              Several Bayesian methods, using the BGLR library:<br>
#'              1.BRR: Bayesian ridge regression: same as rr-blup, but bayesian resolution.<br>
#'              Induces homogeneous shrinkage of all markers effects towards zero with<br>
#'              Gaussian distribution (de los Campos et al, 2013)<br>
#'              2.BL: Bayesian LASSO: uses an exponential prior on marker variances priors,<br>
#'              leading to double exponential distribution of marker effects (Park & Casella 2008)<br>
#'              3.BA: uses a scaled-t prior distribution of marker effects. (Meuwissen et al 2001).<br>
#'              4.BB: Bayes B, uses a mixture of distribution with a point mass at zero<br>
#'              and with a slab of non-zero marker effects with a scaled-t distribution (Habier et al 2011).<br>
#'              5.BC: Bayes C same as Bayes B with a slab with Gaussian distribution.<br>
#'              A more detailed description of these methods can be found in <br>
#'              Perez & de los Campos 2014 (http://genomics.cimmyt.org/BGLR-extdoc.pdf).<br>
#'              Four semi-parametric methods:<br>
#'              1.RKHS: reproductive kernel Hilbert space and multiple kernel MRKHS, using BGLR (Gianola and van Kaam 2008).<br>
#'              Based on genetic distance and a kernel function to regulate the distribution of marker effects.<br>
#'              This methods is claimed to be effective for detecting non additive effects.<br>
#'              2.RF: Random forest regression, using randomForest library (Breiman, 2001, Breiman and Cutler 2013).<br>
#'              This methods uses regression models on tree nodes which are rooted in bootstrapping data.<br>
#'              Supposed to be able to capture interactions between markers<br>
#'              3.SVM: support vector machine, run by e1071 library. For details,<br>
#'              see Chang, Chih-Chung and Lin, Chih-Jen: <br>
#'              LIBSVM: a library for Support Vector Machines http://www.csie.ntu.edu.tw/~cjlin/libsvm<br>
#'              4.LightGBM: Light Gradient Boosting Machine, run by LightGBM library developed by Microsoft.<br>
#'              For details, see Jun Yan, Yuetong Xu, at al: LightGBM: <br>
#'              accelerated genomically designed crop breeding through ensemble learning<br>
#' @param ENalpha used for EN predict model,the value range from 0 to 1.
#' @param nIter Number of iterations for the Bayesian model.
#' @param burnIn Number of burn-in for the Bayesian model.
#' @param thin Number of thins for the Bayesian model.
#' @param SVM_cost Cost of constraints violation for SVM (default: 10).
#' @param kernel the kernel used in training and predicting.
#'               You might consider changing some of the following parameters,
#'               depending on the kernel type.
#'               linear, polynomial, radial basis, sigmoid.
#' @param gamma Parameter needed for all kernels except linear (default: 0.001).
#' @param GBM_params A list of parameters for LightGBM.
#'                   See LightGBM documentation for details.
#' @param GBM_rounds Number of training rounds for LightGBM.
#' @param Envs Assign variable column names of Envs. Default is Envs.
#' @param line_code Assign variable column names of line_code Default is line_Code.
#'
#' @return The class MMPrdM returns a list containing:
#'         The Predicted Values of each line in Multiple Environments.
#'
#' @export
#' @importFrom stats cor lm dist na.omit aggregate median complete.cases
#' @importFrom dplyr %>% mutate group_by row_number summarize_all
#' @importFrom grDevices heat.colors

#' @examples
#' \dontrun{
#' out<-MMPrdM(pheno=pheno,geno=geno,env=env_info,para=envMeanPara)
#' }
MMPrdM<-function(pheno,geno,env,para,Para_Name,model,depend=NULL,reshuffle=NULL,
                 methods=NULL,ENalpha=NULL,SVM_cost=NULL,gamma=NULL,kernel=NULL,
                 Envs=NULL,line_code=NULL,
                 nIter=NULL,burnIn=NULL,thin=NULL,GBM_params=NULL,GBM_rounds=NULL){
  if(is.null(depend)){
    depend="Norm"
  }
  if(is.null(ENalpha)){
    ENalpha=0.5
  }
  if(is.null(methods)){
    methods="RM.G"
  }
  if(is.null(reshuffle)){
    reshuffle= 5
  }
  if(is.null(kernel)){
    kernel="radial"
  }
  if(is.null(Envs)){
    Envs=Envs
  }
  if(is.null(line_code)){
    line_code=line_code
  }
  if(is.null(nIter)){nIter=5000}
  if(is.null(burnIn)){burnIn=1000}
  if(is.null(thin)){thin=3}
  if(is.null(SVM_cost)){SVM_cost=10}
  if(is.null(gamma)){gamma=0.001}
  if(is.null(GBM_rounds)){GBM_rounds=100L}
  removeSaveAt <- function(saveAt) {
    Sys.sleep(2)
    files <- list.files('.', paste0("^", saveAt, '.*(lambda|mu|varE|varB|varU|ScaleBayesA|parBayesB|parBayesC)\\.dat$'))
    tryCatch({
      invisible(lapply(files, file.remove))
    }, error = function(e) {

    })
  }

  enp=which(colnames(para) == Para_Name);
  para.name=colnames(pheno)[-1];

  para=para[para$env_code%in%para.name,];
  env=env[env$env_code%in%para.name,];
  order=match(para$env_code,env$env_code);
  env=env[order,];
  pheno.1=pheno
  keep=dim(pheno.1)[2]-rowSums(is.na(pheno.1))
  pheno=pheno.1
  n.line=dim(pheno)[1];
  n.para=dim(pheno)[2]-1;
  #print(head(pheno))

  coloo=heat.colors(n.para);

  # Use individuals that are present within all environments as test
  # foreach
  id.T <- numeric(0)
  id.V <- numeric(0)
  for (i in 1:n.line) {
    # if NA
    if (any(is.na(pheno[i,-1])) ) {
      id.V <- c(id.V, i)
    }
    else{
      id.T <- c(id.T, i)
    }
  }
  if(depend == "Norm"){
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
      #Assumption 1
      pheno_1 <- pheno[id.T,]
      for(j in 1:length(id.T))
      {

        x1=para[,enp];
        y1=as.vector(t(pheno_1[j,-c(1)]));

        coe=lm(y~x,data=data.frame(x=x1,y=y1));
        inter=summary(coe)$coefficients[1,1]
        slop=summary(coe)$coefficients[2,1];
        intercept=c(intercept,inter);
        slope=c(slope,slop);
      }
      #genotype.match=match(pheno[,1],geno[,1])
      #genotype_1=geno[genotype.match,];
      genotype.impute=rrBLUP::A.mat(geno[,-1],max.missing=0.5,impute.method="mean",return.imputed=T);
      SFI=cbind(geno[,1],genotype.impute$imputed);

      genotype=matrix(suppressWarnings(as.numeric(unlist(SFI))),nrow=nrow(SFI))
      Marker=genotype[,-1];

      intercept.hat=numeric();
      slope.hat=numeric();
      cor.within=matrix(999,reshuffle,n.para);
      cor.all=numeric();
      Prd<-data.frame(matrix(nrow = reshuffle*length(id.V), ncol = n.para))
      colnames(Prd)<-as.vector(para[,2])
      #print(head(Prd))
      for(i in 1:reshuffle){
        Obs<-data.frame()
        yhat.whole.cross=numeric()
        yobs.whole.cross=numeric()
          #id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
          if(model == "rrBLUP"){
            y0<-intercept
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            ans<-rrBLUP::mixed.solve(y=yNa[id.T],Z=genotype[id.T,-1])
            e<-as.matrix(ans$u)
            G.pred<-Marker[id.V,]
            y_pred<-as.matrix(G.pred) %*% e
            GEBV.inter<-c(y_pred[,1])+c(ans$beta)

            y0<-slope
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            ans<-rrBLUP::mixed.solve(y=yNa[id.T],Z=genotype[id.T,-1])
            e<-as.matrix(ans$u)
            G.pred<-Marker[id.V,]
            y_pred<-as.matrix(G.pred) %*% e

            GEBV.slope<-c(y_pred[,1])+c(ans$beta);
            #print("Model")
          }
        else if(model == "LASSO"){
            y0<-intercept
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            cv.fit<-glmnet::cv.glmnet(y=yNa[id.T],x=genotype[id.T,-1],
                                      alpha=1)
            lambda_min <- cv.fit$lambda.min
            G.pred<-Marker[id.V,]
            GEBV.inter<- stats::predict(cv.fit,newx=G.pred,s=c(lambda_min))
            #coef<-coef(cv.fit)
            #e=as.matrix(coef@x[-1])
            # G.pred=Marker[id.V,]
            #selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            # G.pred_LASSO <- G.pred[, selected_features]
            # y_pred=as.matrix(G.pred_LASSO) %*% e
            #GEBV.inter=c(y_pred[,1])+c(coef@x[1]);

            y0<-slope
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            cv.fit<-glmnet::cv.glmnet(y=yNa[id.T],x=genotype[id.T,-1],
                                      alpha=1)
            lambda_min <- cv.fit$lambda.min
            G.pred<-Marker[id.V,]
            GEBV.slope<- stats::predict(cv.fit,newx=G.pred,s=c(lambda_min))
            #coef<-coef(cv.fit)
            #e=as.matrix(coef@x[-1])
            #G.pred=Marker[id.V,]
            #selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            #G.pred_LASSO <- G.pred[, selected_features]
            #y_pred=as.matrix(G.pred_LASSO) %*% e
            #GEBV.slope=c(y_pred[,1])+c(coef@x[1]);
          }else if(model == "RR"){
            y0<-intercept
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            cv.fit<-glmnet::cv.glmnet(y=yNa[id.T],x=genotype[id.T,-1],
                                      alpha=0)
            lambda_min <- cv.fit$lambda.min
            G.pred<-Marker[id.V,]
            GEBV.inter<- stats::predict(cv.fit,newx=G.pred,s=c(lambda_min))
            #coef<-coef(cv.fit)
            #e=as.matrix(coef@x[-1])
            # G.pred=Marker[id.V,]
            #selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            # G.pred_RR <- G.pred[, selected_features]
            # y_pred=as.matrix(G.pred_RR) %*% e
            #GEBV.inter=c(y_pred[,1])+c(coef@x[1]);

            y0<-slope
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            cv.fit<-glmnet::cv.glmnet(y=yNa[id.T],x=genotype[id.T,-1],
                                      alpha=0)
            lambda_min <- cv.fit$lambda.min
            G.pred<-Marker[id.V,]
            GEBV.slope<- stats::predict(cv.fit,newx=G.pred,s=c(lambda_min))
            #coef<-coef(cv.fit)
            #e=as.matrix(coef@x[-1])
            #G.pred=Marker[id.V,]
            #selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            #G.pred_RR <- G.pred[, selected_features]
            #y_pred=as.matrix(G.pred_RR) %*% e
            #GEBV.slope=c(y_pred[,1])+c(coef@x[1]);
          }
          else if(model == "EN"){
            y0<-intercept
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            cv.fit<-glmnet::cv.glmnet(y=yNa[id.T],x=genotype[id.T,-1],
                                      alpha=ENalpha)
            lambda_min <- cv.fit$lambda.min
            G.pred<-Marker[id.V,]
            GEBV.inter<- stats::predict(cv.fit,newx=G.pred,s=c(lambda_min))
            #coef<-coef(cv.fit)
            #e=as.matrix(coef@x[-1])
            # G.pred=Marker[id.V,]
            #selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            # G.pred_EN <- G.pred[, selected_features]
            # y_pred=as.matrix(G.pred_EN) %*% e
            #GEBV.inter=c(y_pred[,1])+c(coef@x[1]);

            y0<-slope
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            cv.fit<-glmnet::cv.glmnet(y=yNa[id.T],x=genotype[id.T,-1],
                                      alpha=ENalpha)
            lambda_min <- cv.fit$lambda.min
            G.pred<-Marker[id.V,]
            GEBV.slope<- stats::predict(cv.fit,newx=G.pred,s=c(lambda_min))
            #coef<-coef(cv.fit)
            #e=as.matrix(coef@x[-1])
            #G.pred=Marker[id.V,]
            #selected_features <- which(as.vector(coef(cv.fit)[-1, ]) != 0)  # 获取非零系数的特征索引
            #G.pred_EN <- G.pred[, selected_features]
            #y_pred=as.matrix(G.pred_EN) %*% e
            #GEBV.slope=c(y_pred[,1])+c(coef@x[1]);
          } else if(model == "RF"){
            y0<- intercept
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            GEBV.inter<-randomForest::randomForest(genotype[id.T,-1], y=yNa[id.T],
                                                   xtest=Marker[id.V,])$test$predicted

            y0<- slope
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            GEBV.slope<-randomForest::randomForest(genotype[id.T,-1], y=yNa[id.T],
                                                   xtest=Marker[id.V,])$test$predicted
          }
          else if(model == "SVM"){
            y0<- intercept
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            model.inter<-e1071::svm(genotype[id.T,-1], y=yNa[id.T],
                                    method="nu-regression",
                                    kernel="linear",cost=SVM_cost,
                                    gamma=gamma,max_iter=100)
            GEBV.inter<-stats::predict(model.inter,Marker[id.V,])

            y0<- slope
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            model.slope<-e1071::svm(genotype[id.T,-1], y=yNa[id.T]*10^9,
                                    method="nu-regression",
                                    kernel="linear",cost=SVM_cost,
                                    gamma=gamma,max_iter=100)
            GEBV.slope<-stats::predict(model.slope,Marker[id.V,])/10^9
          }
          else if(model == "BA" | model == "BC" | model == "BL" | model == "BRR"
                  | model == "BB" ){
            y0=intercept;
            GENO = genotype[,-1]
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
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
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
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
            y0<-intercept
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            GENO<-as.data.frame(geno[,-1])
            rownames(GENO)<-geno[,1]
            A <- rrBLUP::A.mat(as.data.frame(GENO))
            #rownames(A)<-geno[,1]
            dataF <- data.frame(genoID = rownames(GENO), intcp = yNa)
            #A<-rrBLUP::A.mat(as.data.frame(geno[id.V,-1]))
            #dataF=data.frame(genoID=rownames(geno[id.V,-1]),intcp=y0[id.V])
            MODEL<-rrBLUP::kin.blup(data=dataF,geno="genoID",pheno="intcp", GAUSS=FALSE, K=A,
                                    PEV=TRUE,n.core=1,theta.seq=NULL)
            GEBV.inter<-MODEL$pred[id.V]
            #GEBV.inter<-MODEL$pred


            y0<-slope
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            GENO<-as.data.frame(genotype[,-1])
            rownames(GENO)<-geno[,1]
            A <- rrBLUP::A.mat(as.data.frame(GENO))
            rownames(A)<-geno[,1]
            dataF <- data.frame(genoID = rownames(GENO), intcp = yNa)
            #A<-rrBLUP::A.mat(as.data.frame(geno[id.V,-1]))
            #dataF=data.frame(genoID=rownames(geno[id.V,-1]),intcp=y0[id.V])
            MODEL<-rrBLUP::kin.blup(data=dataF,geno="genoID",pheno="intcp", GAUSS=FALSE, K=A,
                                    PEV=TRUE,n.core=1,theta.seq=NULL)
            GEBV.slope<-MODEL$pred[id.V]
            #GEBV.slope<-MODEL$pred

          }
          else if(model =="RKHS" | model =="MKRKHS"){
            y0=intercept;
            GENO = genotype[,-1]
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA

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
            removeSaveAt(saveAt)

            y0=slope;
            GENO = genotype[,-1]
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA

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
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            dtrain <- lightgbm::lgb.Dataset(genotype[id.T,-1], label = yNa[id.T])
            dtest <- lightgbm::lgb.Dataset.create.valid(dtrain, genotype[id.V,-1], label = yNa[id.V])
            valids <- list(test = dtest)
            model.inter <- lightgbm::lgb.train(params = params,data = dtrain
                                               ,valids = valids,nrounds = GBM_rounds,verbose = -1)
            GEBV.inter <- stats::predict(model.inter, Marker[id.V,])
            #print(head(GEBV.inter))

            y0<-slope
            yNa<- c()
            yNa[id.T]<-y0
            yNa[id.V]<-NA
            dtrain <- lightgbm::lgb.Dataset(genotype[id.T,-1], label = yNa[id.T])
            dtest <- lightgbm::lgb.Dataset.create.valid(dtrain, genotype[id.V,-1], label = yNa[id.V])
            valids <- list(test = dtest)
            model.slope <- lightgbm::lgb.train(params = params,data = dtrain,
                                               valids = valids,nrounds = GBM_rounds,verbose = -1)
            GEBV.slope <- stats::predict(model.slope, Marker[id.V,])
            #print(head(GEBV.slope))
            #rm(ls(model.slope,model.slope))
            #return(model.inter)
          }
          yhat.para=matrix(999,length(id.V),n.para)
          yobs.para=matrix(999,length(id.V),n.para)
          for(j in 1:n.para)
          {
            yhat<-GEBV.inter+GEBV.slope*para[j,enp];
            yobs<-pheno[id.V,j+1];
            #return(list(yhat,yhat.para[,1]))
            yhat.para[,j]<-yhat
            yobs.para[,j]<-yobs
          }
          yhat.whole.cross<-rbind(yhat.whole.cross,yhat.para)
          yobs.whole.cross<-rbind(yobs.whole.cross,yobs.para)
          s<-length(id.V)*(i-1)+1
          t<-length(id.V)*i
          Prd[c(s:t),]<-yhat.whole.cross
      }
      Prd$line_code<-rep(pheno[id.V,1],reshuffle)
      Prd_mean <- Prd %>% group_by(line_code) %>% summarize_all(mean)
      Prd_mean <- as.data.frame(Prd_mean)
      #Correlation within environment 50 times
      return(Prd=Prd_mean)
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

          yhat.whole=numeric();
          yobs.whole=numeric();

            #id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];

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
                                      kernel="radial",cost=SVM_cost,gamma=gamma)
              GEBV.inter<-stats::predict(model.inter,Marker[id.V,])

              y0=slope;
              model.slope<-e1071::svm(genotype[id.T,-1], y=y0[id.T], method="nu-regression",
                                      kernel="radial",cost=SVM_cost,gamma=gamma)
              GEBV.slope<-stats::predict(model.slope,Marker[id.V,])
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
              MODEL=rrBLUP::kin.blup(data=dataF,geno="genoID",pheno="intcp",
                                     GAUSS=FALSE, K=A,
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
  }else if(depend =="PEI"){
    genotype.match=match(pheno[,1],geno[,1])
    genotype_1=geno[genotype.match,];
    genotype.impute=rrBLUP::A.mat(genotype_1[,-1],max.missing=0.5,
                                  impute.method="mean",return.imputed=T);
    SFI=cbind(genotype_1[,1],genotype.impute$imputed);
    genotype=matrix(suppressWarnings(as.numeric(unlist(SFI))),nrow=nrow(SFI))
    Marker=genotype[,-1];
    n.marker=dim(Marker)[2];
    if(methods=="RM.G"){
      cor.within<-matrix(999,reshuffle,n.para)
      cor.all<-numeric()
      if(model == "rrBLUPJ" | model == "RRJ"){
        for(i in 1:reshuffle)
        {
          yhat.whole.cross=numeric();
          yobs.whole.cross=numeric();


            #id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
            ##marker effect###
            effect=matrix(999,n.marker,n.para);
            intercept=matrix(999,1,n.para)

            for(k in 1:n.para)
            {
              id.T <- c(1:n.line)[!is.na(pheno[,1+k])]
              id.V <- c(1:n.line)[!is.na(pheno[,1+k])]
              if(model == "rrBLUPJ"){
                fit=rrBLUP::mixed.solve(pheno[id.T,1+k],Z=Marker[id.T,])
                effect[,k]=fit$u
                intercept[,k]=fit$beta
              }
              if(model == "RRJ"){
                cv.fit<-glmnet::cv.glmnet(y=pheno[id.T,1+k],x=Marker[id.T,],
                                          alpha=0)
                lambda_min <- cv.fit$lambda.min
                cv.fit<-glmnet::cv.glmnet(y=pheno[id.T,1+k],x=Marker[id.T,],
                                          alpha=0,lambda=lambda_min)
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
        PrdF<-data.frame()
        for(r in 1:reshuffle){
          yhat.whole.cross<-data.frame()
          yobs.whole.cross<-data.frame()

          id.V<-c()
          id.T<-c()

          # Use individuals that are present within all environments as test
          # foreach
          for (i in 1:n.line) {
            # if NA
            if (any(is.na(pheno[i,-1])) ) {
              id.V <- c(id.V, i)
            }
            else{
              id.T <- c(id.T, i)
            }
          }
          print(length(id.V))
          if(length(id.V) != dim(pheno)[1]){
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
            #pheno3
            y <- pheno
            ys <- y[order(factor(y$line_code, levels = geno[,1])), ]
            yl <- ys %>%tidyr::pivot_longer(cols = -line_code,
                                            names_to = "Envs",
                                            values_to = "Obs")
            yl<-as.data.frame(yl)
            y<-yl[order(yl[,2]),]

            #pheno1
            yNa  <- pheno[id.T,]
            yNas <- yNa[order(factor(yNa$line_code, levels = geno[id.T,1])), ]

            yNal <- yNas %>%tidyr::pivot_longer(cols = -line_code,
                                                names_to = "Envs",
                                                values_to = "Obs")
            yNal<-as.data.frame(yNal)
            yNa<-yNal[order(yNal[,2]), ]
          }
          else if(length(id.V) == dim(pheno)[1]){
            #GENO-ALL
            n <- nrow(as.matrix(geno[,-1]))
            m <- ncol(as.matrix(geno[,-1]))
            n.envs_all<- dim(pheno)[2]-1
            GENO <- matrix(0, nrow = n * n.envs_all, ncol = m * n.envs_all)
            for (i in 1:n.envs_all) {
              row_start <- (i - 1) * n + 1
              row_end <- i * n
              col_start <- (i - 1) * m + 1
              col_end <- i * m
              GENO[row_start:row_end, col_start:col_end] <- as.matrix(geno[,-1])
            }
            #PHENO-ALL
            y <- pheno
            ys <- y[order(factor(y$line_code, levels = geno[,1])), ]
            yl <- ys %>%tidyr::pivot_longer(cols = -line_code,
                                            names_to = "Envs",
                                            values_to = "Obs")
            yl<-as.data.frame(yl)
            y<-yl[order(yl[,2]),]

            #PHENO for Train
            yNa<-na.omit(y)

            #GENO for Train
            M <- matrix(0, nrow = n * (n.envs_all-1), ncol = m * n.envs_all)
            for (i in 1:(n.envs_all-1)) {
              row_start <- (i - 1) * n + 1
              row_end <- i * n
              col_start <- (i - 1) * m + 1
              col_end <- i * m
              M[row_start:row_end, col_start:col_end] <- as.matrix(geno[,-1])
            }
            id <-colnames(pheno)[which(is.na(pheno), arr.ind=TRUE)[1,2]]
          }
            if (model == "RF"){
              Prd<-randomForest::randomForest(M, y=yNa[,ncol(yNa)],xtest=GENO)$test$predicted
              PrdM<-data.frame(y, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = "line_code",
                                                      names_from ="Envs",
                                                      values_from = "Prd")
              if((length(id.V) != dim(pheno)[1])){
                PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% pheno[id.V,1],]
              }
              else if(length(id.V) == dim(pheno)[1]){
                PrdM_wide<- data.frame(line_code=PrdM_wide$line_code,Prd=PrdM_wide[[id]])
              }
              PrdM_wide<-as.data.frame(PrdM_wide)
            }
            if (model == "rrBLUP"){
              ans<-rrBLUP::mixed.solve(y=yNa[,ncol(yNa)],Z=M)
              e<-as.matrix(ans$u)
              G.pred<-GENO
              y_pred<-as.matrix(G.pred) %*% e
              Prd<-c(y_pred[,1])+c(ans$beta)
              PrdM<-data.frame(y, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(
                id_cols = "line_code",
                names_from ="Envs",
                values_from = "Prd")
              PrdM_wide<-as.data.frame(PrdM_wide)
              if((length(id.V) != dim(pheno)[1])){
                PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% pheno[id.V,1],]
              }
              else if(length(id.V) == dim(pheno)[1]){
                PrdM_wide<- data.frame(line_code=PrdM_wide$line_code,Prd=PrdM_wide[[id]])
              }
              PrdM_wide<-as.data.frame(PrdM_wide)
              return(PrdM_wide)
            }
            else if (model == "RR"){
              cv.fit <- glmnet::cv.glmnet(M,yNa[,ncol(yNa)],family="gaussian",alpha=0)
              lambda_min <- cv.fit$lambda.min
              Prd<- stats::predict(cv.fit,newx=GENO,s=c(lambda_min))
              PrdM<-data.frame(y, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(
                id_cols = "line_code",
                names_from ="Envs",
                values_from = "Prd")
              if((length(id.V) != dim(pheno)[1])){
                PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% pheno[id.V,1],]
              }
              else if(length(id.V) == dim(pheno)[1]){
                PrdM_wide<- data.frame(line_code=PrdM_wide$line_code,Prd=PrdM_wide[[id]])
              }
              PrdM_wide<-as.data.frame(PrdM_wide)
            }
            else if (model == "LASSO"){
              cv.fit <- glmnet::cv.glmnet(M,yNa[,ncol(yNa)],family="gaussian",alpha=1)
              lambda_min <- cv.fit$lambda.min
              Prd<- stats::predict(cv.fit,newx=GENO,s=c(lambda_min))
              PrdM<-data.frame(y, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(
                id_cols = "line_code",
                names_from ="Envs",
                values_from = "Prd")
              if((length(id.V) != dim(pheno)[1])){
                PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% pheno[id.V,1],]
              }
              else if(length(id.V) == dim(pheno)[1]){
                PrdM_wide<- data.frame(line_code=PrdM_wide$line_code,Prd=PrdM_wide[[id]])
              }
              PrdM_wide<-as.data.frame(PrdM_wide)
            }
            else if (model == "EN"){
              cv.fit <- glmnet::cv.glmnet(M,yNa[,ncol(yNa)],family="gaussian",alpha=ENalpha)
              lambda_min <- cv.fit$lambda.min
              Prd<- stats::predict(cv.fit,newx=GENO,s=c(lambda_min))
              PrdM<-data.frame(y, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(
                id_cols = "line_code",
                names_from ="Envs",
                values_from = "Prd")
              if((length(id.V) != dim(pheno)[1])){
                PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% pheno[id.V,1],]
              }
              else if(length(id.V) == dim(pheno)[1]){
                PrdM_wide<- data.frame(line_code=PrdM_wide$line_code,Prd=PrdM_wide[[id]])
              }
              PrdM_wide<-as.data.frame(PrdM_wide)
            }
            else if (model == "GBLUP"){
              GENO<-as.data.frame(GENO)
              rowname<-paste(y[,1],y[,2], sep="-")
              rownames(GENO)<-rowname
              A<-rrBLUP::A.mat(GENO)
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
                id_cols = "line_code",
                names_from ="Envs",
                values_from = "Prd")
              PrdM_wide<-as.data.frame(PrdM_wide)
              if((length(id.V) != dim(pheno)[1])){
                PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% pheno[id.V,1],]
              }
              else{(length(id.V) != dim(pheno)[1])
                PrdM_wide<- data.frame(line_code=PrdM_wide$line_code,Prd=PrdM_wide[[id]])
              }
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
                id_cols = "line_code",
                names_from ="Envs",
                values_from = "Prd")
              PrdM_wide<-as.data.frame(PrdM_wide)
              if((length(id.V) != dim(pheno)[1])){
                PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% pheno[id.V,1],]
              }
              else if(length(id.V) == dim(pheno)[1]){
                PrdM_wide<- data.frame(line_code=PrdM_wide$line_code,Prd=PrdM_wide[[id]])
              }
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
                id_cols = "line_code",
                names_from ="Envs",
                values_from = "Prd")
              if((length(id.V) != dim(pheno)[1])){
                PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% pheno[id.V,1],]
              }
              else if(length(id.V) == dim(pheno)[1]){
                PrdM_wide<- data.frame(line_code=PrdM_wide$line_code,Prd=PrdM_wide[[id]])
              }
              removeSaveAt(saveAt)
            }
            else if (model == "SVM"){
              model_SVM <- e1071::svm(M,yNa[,ncol(yNa)],method="nu-regression",
                                      kernel=kernel,cost=SVM_cost,gamma=gamma)
              Prd <- stats::predict(model_SVM,GENO)
              PrdM<-data.frame(y, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = "line_code",
                                                      names_from ="Envs",
                                                      values_from = "Prd")
              if((length(id.V) != dim(pheno)[1])){
                PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% pheno[id.V,1],]
              }
              else if(length(id.V) == dim(pheno)[1]){
                PrdM_wide<- data.frame(line_code=PrdM_wide$line_code,Prd=PrdM_wide[[id]])
              }
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
              dtest <- lightgbm::lgb.Dataset.create.valid(dtrain, GENO, label = y[,ncol(y)])
              valids <- list(test = dtest)
              model_GBM <- lightgbm::lgb.train(params = params,data = dtrain,
                                               valids = valids,nrounds = GBM_rounds,verbose = -1)
              Prd <- stats::predict(model_GBM, GENO)
              PrdM<-data.frame(y, Prd = as.numeric(Prd))
              PrdM_wide<- PrdM %>% tidyr::pivot_wider(id_cols = "line_code",
                                                      names_from ="Envs",
                                                      values_from = "Prd")
              if((length(id.V) != dim(pheno)[1])){
                PrdM_wide<-PrdM_wide[PrdM_wide$line_code %in% pheno[id.V,1],]
              }
              else if(length(id.V) == dim(pheno)[1]){
                PrdM_wide<- data.frame(line_code=PrdM_wide$line_code,Prd=PrdM_wide[[id]])
              }
              PrdM_wide<-as.data.frame(PrdM_wide)
            }
            PrdF<-rbind(PrdF,PrdM_wide)
            }
        Prd_mean<-PrdF %>% group_by(line_code) %>% summarize_all(mean)
        Prd_mean<-as.data.frame(Prd_mean)
        return(Prd_mean)
        #Correlation within environment 50 times
      }
      if(methods=="RM.GE"){
        cor.within=matrix(999,reshuffle,n.para);cor.all=numeric();
        for(i in 1:reshuffle)
        {
          #cross = sample(rep(1:fold,each=ceiling(n.line/fold)),n.line)
          obs_matrix=numeric();pre_matrix=numeric();

            #id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
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
