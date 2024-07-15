
# MMGS
A comprehensive R package, multiple-environments multiple methods genomic selection (MMGS), developed by Mingjia Zhu, integrates the polygenic environmental interaction (PEI) and Reaction Norm (RE) methods along with 15 prediction models that include difference prediction estimated methods contains parametic, semi-parametric and non-parametric.

## Description
RE model includes four steps: (1) Using CERIS algorithm (Guo, 2021) to identify an environmental index that explained the largest proportion of phenotypic variation. (2) Regressing the observed phenotypes on the identified environmental index to obtain an intercept and a slope estimate for each tested genotype. (3) Treating intercept and slope as new "traits" and perform genomic prediction through ridge regression to predict the intercept and slope for each untested genotype. (4) Predict the phenotypes of the untested genotypes using the predicted intercept and slope and the environmental index value of each environment. Consistent with the RE model, the PEI model starts with identifying key environmental index that best captures the phenotypic variation

Total of these predicted statistical models are classified into three major categories: parametric, semi-parametric, and non-parametric (Admas et al, 2024). The parametric statistical models include mixed linear models like genomic best linear unbiased prediction (G-BLUP) (Vanraden, 2008), BayesA (BA) and BayesB (BB) (Meuwissen et al., 2001), BayesC (BC) (George and McCulloch, 1993), Bayesian ridge regression (BRR) (Erbe et al., 2012), and Bayesian LASSO (BL) (Park and Casella, 2008), least absolute shrinkage and selection operator (LASSO) (Usai et al., 2009), ridge regression (RR) (Whittaker et al., 2000), ridge regression best linear unbiased prediction (RR-BLUP) (Meuwissen et al., 2001), and elastic net (EN) (Zou and Hastie, 2005). The semi-parametric method includes the reproducing kernel Hilbert space (RKHS) model and multiple kernel RKHS (MKRKHS) (Gianola et al., 2006). The non-parametric method comprises support vector machine (SVM) (Maenhout et al., 2007), and random forest (RF) (Chen and Ishwaran, 2012), and gradient boosting machine (GBM) (Li et al., 2018). 
## Installation

You can install the package from CRAN using the following command:
From Github:
```R
devtools::install_github("Ryougi-yukiro/MMGS")
```

## Usage
Provide examples of how to use your package. Include code snippets and brief explanations to demonstrate the key features and functionalities.
You can also provide links to additional resources or documentation.

### Step.1 Load packages and data
We have built-in data from a hybrid population that includes environmental data from multiple locations, filtered genotype data, and flowering-related phenotypic data. This dataset is smaller and easier for beginners to understand how the package is used.

```R
#Load the required packages
library("MMGS")
library("dplyr")# used for data reshape and melt

data(trait)
data(geno)
data(env_info)
data("PTT_PTR")
```


### Step.2 data analysis 
This step explores the basic attributes of the data situation and provides pre-processing for subsequent analysis.
```R
env_trait<-env_trait_calculate(data=trait,trait="FTgdd",env="env_code")
LbyE<-LbyE_calculate(data=trait,trait="FTgdd",env="env_code",line="line_code")
LbyE_corrplot(LbyE=LbyE)
etl<-LbyE_Reshape(data=env_trait,env="env_code",LbyE=LbyE)
etl_plotter(data=etl,trait=env_trait)
Regression<-Reg(LbyE = LbyE, env_trait = env_trait)
#Reg_plotter(Reg = Regression)
result<-line_trait_mean(data=trait,trait="FTgdd",mean=env_trait,LbyE=LbyE,row=2)
MSE<-result[[1]]
ltm<-result[[2]]
#mse_plotter(MSE)
#Mean_trait_plot(Regression,MSE)
```

### Step.3 Search env index
This step aims to find the most relevant environmental factors to provide a solid basis for subsequent predictions
```R
Paras <- colnames(PTT_PTR)[-c(1:4)]
#windows-search
pop_cor<-Exhaustive_search(data=env_trait, env_paras=PTT_PTR, searching_daps=122,
                           p=1, dap_x=122,dap_y=122,LOO=0,Paras=Paras)
#plot
#Exhaustive_plotter(Correlation=pop_cor,dap_x=122, dap_y=122,p=1,Paras=Paras)

#correlation
envMeanPara<-envMeanPara(data=env_trait, env_paras=PTT_PTR, maxR_dap1=18,
                         maxR_dap2=43, Paras=Paras)
#plot
#envMeanPara_plotter(data=envMeanPara,Paras=Paras)
```

### Step.3 CV
Users can customize the model they need, the function uses the  by default, the given environment parameters can be obtained from the previous results , fold number represents the number of folds, reshuffle represents the number of repetitions. 
```R
#Check pheno
pheno<-LbyE[which(as.character(LbyE$line_code)%in%c("line_code",as.character(geno$line_code))),];
#CV 
out<-MMGP(pheno=pheno, geno=geno, env=env_info,para=envMeanPara, Para_Name=Para[1], depend="PEI",model="BB", kernel="linear", fold=2, reshuffle=5, methods="RM.G")
#result
#> mean(out[[3]])
#[1] 0.8728506
#> apply(out[[2]],2,mean)
#     PR12      IA14      PR11      IA13     PR14S      KS11      KS12 
#0.5418663 0.3868576 0.5381628 0.4759335 0.4871427 0.6219213 0.6380658
#head(out[[1]])
#       obs      pre     col para
#1 1595.988 1588.782 #FF0000 PR12
#2 1512.918 1576.437 #FF0000 PR12
```
Correlation here refers to the correlation between the predicted phenotypes and the actual phenotypes of the environment, not the breeding values, so please do your own calculations first if needed (before the R package is updated).

### Others function Example
```R
pheno<-LbyE
pheno$PR11 <-NA
#linear radial polynomial linear
#library(dplyr)
for( i in envMeanPara$env_code){
  pheno<-LbyE
  pheno[["KS12"]]<-NA
out<-MMPrdM(pheno=pheno, geno=geno,env=env_info,para=envMeanPara,
            Para_Name=c("PTS"), depend="PEI",
            SVM_cost = 1,gamma=10,kernel="linear",fixed=T,
            model="SVM",reshuffle=1,methods="RM.G")
(cor<-cor(out[,2],LbyE[["KS12"]]))
print(paste(i," : ",cor))
}
```

## Documentation
See full documentation from original repository

## Command line interface
* env_trait_calculate
* envMeanPara
* envMeanPara_plotter
* etl_calculate
* etl_plotter
* Exhaustive_plotter
* Exhaustive_search
* MMGS
* h2_rrBLUP
* LbyE_calculate
* LbyE_corrplot
* line_trait_mean
* ltm_plotter
* Mean_trait_plot
* mse_plotter
* prdM_plotter
* Reg
* Reg_plotter
* Slope_Intercept

## Implementation notes
MMGS is a collection of tools for cross-environmental genome-wide selection prediction that integrates most genome-wide prediction models, both parametric and non-parametric. You can input your own collected data against sample data and get the results you want directly through the built-in functions of the toolkit, which requires no additional statistical knowledge or coding skills and is somewhat user-friendly because it saves users from having to search for various tools and apply them to cross-environmental prediction.

## License
[GPL-3](https://www.gnu.org/licenses/quick-guide-gplv3.html)

## Reference
[Li X, Guo T, Wang J, et al. An integrated framework reinstating the environmental dimension for GWAS and genomic selection in crops[J]. Molecular Plant](https://www.sciencedirect.com/science/article/pii/S167420522100085X)

[Jarquín D, Crossa J, Lacaze X, et al. A reaction norm model for genomic selection using high-dimensional genomic and environmental data[J]. Theoretical and applied genetics, 2014, 127: 595-607.](https://link.springer.com/article/10.1007/s00122-013-2243-1)








