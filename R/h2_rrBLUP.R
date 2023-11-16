#' Title h2_rrBLUP caculate
#'
#' @param trait pheno_trait
#' @param geno geno_type_txt,0,1,2type
<<<<<<< HEAD
#' @param envs envs you want to calculate
=======
#' @param envs
>>>>>>> 23b0efae316a77181b94cb5eff39d6257393dcf5
#'
#' @return h2 in different envs
#' @export
#'
<<<<<<< HEAD
#' @examples \dontrun{h2<-h2_rrBLUP(trait=trait,geno=geno,envs="env_code")}
=======
#' @examples h2<-h2_rrBLUP(trait=trait,geno=geno,envs="env_code")
>>>>>>> 23b0efae316a77181b94cb5eff39d6257393dcf5
h2_rrBLUP<-function(trait,geno,envs){
  c.z.hglm <- function(kin){
    relmat <- kin
    relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]
    svd <- svd(relmat)
    Z <- svd$u %*% diag(sqrt(svd$d))
    return(Z)
  }
  h2_rrblup <- function(y,Z){
    est <- rrBLUP::mixed.solve(y=y,Z = Z,X = as.matrix(rep(1,length(y))),
                               method="reml",bounds=c(1e-09,1e+09),SE=FALSE)
    h2<-est$Vu/(est$Vu+est$Ve)
    return(h2)
  }
  kin_ibs_matrix2<-as.matrix(rrBLUP::A.mat(geno[,-1]))
  Z <- c.z.hglm(kin_ibs_matrix2)
  d1<-data.frame()
  for(i in unique(trait[[envs]])){
    data<-subset(trait,trait[[envs]] == i )
<<<<<<< HEAD
    y=as.numeric(data[[envs]])
=======
    y=as.numeric(data[,5])
>>>>>>> 23b0efae316a77181b94cb5eff39d6257393dcf5
    a<-h2_rrblup(y,Z)
    a<-c(i,a)
    d1<-rbind(d1,a)
  }
  colnames(d1)<-c("envs","h2")
  return(d1)
}
