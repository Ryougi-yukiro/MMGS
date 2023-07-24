# ggGSE
The R package ggGSE is developed by . You can use this packages for Genomic selection across different envs

## Description

A  devlopment R packages used for GS within diff envs. 
CV based on norm reaction  have been modified in our R packages.
Almost of Genomic prediction has been accommodated.
Others like Deep_GS, due to its indv docker, we chose to drop the introduction of this model in order to ensure the generalizability and user-friendliness of this R package;
For DNNGP, One, because the core code is compiled and his CNN network model is not available, and two, because his input file is a PCA file.

## Installation

You can install the package from CRAN using the following command:
From Github:
```R
devtools::install_github("Ryougi-yukiro/ggGSE")
```

## Usage
Provide examples of how to use your package. Include code snippets and brief explanations to demonstrate the key features and functionalities.
You can also provide links to additional resources or documentation.

### Step.1 Load packages and data
```R
#Load the required packages
library("ggGSE")
library("dplyr")# used for data reshape and melt

#Load the input data
trait<-read.table(file="Trait_records.txt",header=T)
env_info<-read.table(file="Env_meta_table.txt",header=T)
PTT_PTR <- read.table(file="7Envs_envParas_DAP122.txt", header = T , sep = "\t")
geno <- read.table(file="Genotype.txt", header = T , sep = "\t")
```
### Step.2 estimate h2 within different Envs
```R
h2<-h2_rrBLUP(trait=trait,geno=geno,envs="env_code")
#if you need estiamte H2 plz use asreml or VCA package
#example code as follow:
trait$Loc <- as.factor(trait$env_code)
trait$GDD <- as.numeric(as.character(trait$PH))
trait$Entry <- gsub("E", "", trait$line_code)

Loc=as.character(trait$Loc);
Entry=as.character(trait$Entry);
Rep=as.character(trait$pop_code);
y=trait$GDD;

mydataframe=data.frame(y,Loc,Entry,Rep);
fit_r <-remlMM(y ~ (Loc)+(Entry)+(Loc/Rep)+(Loc*Entry), mydataframe)
GxE_Var=fit_r$aov.tab[5,2];
G_Var=fit_r$aov.tab[3,2];
Error_Var=fit_r$aov.tab[6,2];
fit_a <-anovaMM(y ~ Loc+Loc/Rep+Entry+Loc*Entry, mydataframe)
fit_a
EN=unique(Loc);
GenoVarInd=numeric();
for(i in 1:7)
{
  data1=mydataframe[mydataframe$Loc==EN[1],];
  data1=data1[which(!(is.na(data1$y))),];
  fit1 <-remlMM(y ~ (Entry)+(Rep), data1)
  print(fit1)
  #Variance component for genotypes###
  GenoVar=fit1$aov.tab[2,2];
  GenoVarInd=c(GenoVarInd,GenoVar)
}
#V
V=sum((sqrt(GenoVarInd)-mean(sqrt(GenoVarInd)))^2)/(length(EN)-1)
V.percent=V/GxE_Var;
V
V.percent
#L
L=GxE_Var-V;
L.percent=L/GxE_Var;
L
L.percent

#Line/Entry mean heritability
h2=G_Var/(G_Var+GxE_Var/length(EN)+Error_Var/(length(EN)*2));
h2
#Pooled genetic correlation
rg=G_Var/(G_Var+L)
rg
```
### Step.2 data analysis 
```R
env_trait<-env_trait_calculate(data=trait,trait="FTgdd",env="env_code")
LbyE<-LbyE_calculate(data=trait,trait="FTgdd",env="env_code",line="line_code")
#envs similarity, plot by ggplot2 like heatmap
LbyE_corrplot(LbyE=LbyE,cor_type="phetamp",color=c("blue","white","red"))

#Get lines mean and Q25, Q75 within each envs
etl<-etl_calculate(data=env_info,trait="FTgdd",env="env_code",bycol="lat")
etl_plotter(data=etl,mean=env_trait) # you could design your plotter use more paras

#Get Best Days window and Env Paras
Paras <- c('DL', 'GDD', 'PTT', 'PTR', 'PTS');
#p dap_x dap_y searching_daps according to your data
pop_cor<-Exhaustive_search(data=env_trait, env_paras=PTT_PTR, searching_daps=80,
                           p=1, dap_x=80,dap_y=80,LOO=0,Paras=Paras)
#Result Visualization
Exhaustive_plotter(data=pop_cor,dap_x=80, dap_y=80,p=1)


dap1<-pop_cor[,2][which(pop_cor[,14]==max(abs(pop_cor[,14])))]
dap2<-pop_cor[,3][which(pop_cor[,14]==max(abs(pop_cor[,14])))]

#Get Envs mean within env paras
envMeanPara<-envMeanPara(data=env_trait, env_paras=PTT_PTR, maxR_dap1=18,maxR_dap2=43, Paras=Paras)
```

### Step.3 CV
Users can customize the model they need, the function uses the norm reaction by default, the given environment parameters can be obtained from the previous results , fold number represents the number of folds, reshuffle represents the number of repetitions. In RM.G mode, the available models are rrBLUP,GLUP,LASSO,EN,RR,BA,BB,BC,BL,BRR,RKHS,MKRKHS,SVM,RF and LightGBM.
```R
#Check pheno
pheno<-LbyE[which(as.character(LbyE$line_code)%in%c("line_code",as.character(geno$line_code))),];
#CV 
out<-GE_CV(pheno=pheno, geno=geno, env=env_info,
             para=envMeanPara, Para_Name="PTT", depend="norm",
             model="rrBLUP", fold=2, reshuffle=5, methods="RM.G")
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

### For Deep Learning Model LightGBM

For deep learning models, different choices of hyperparameters can greatly affect the results of model predictions.
Therefore, for samples with different population sizes you need to find different hyperparameters to determine the fitted results.

Although we write the LightGBM model into the R package with a pre-set hyperparameter, this is only a not-so-bad parameter from different datasets and does not represent the best prediction result of the model, so when using this model, make sure to customize your hyperparameter to ensure the best prediction result

Below is our pre-set list of hyperparameters, please refer to the official LightGBM documentation for the meaning of the parameters.

```R
            params <- list(boosting="gbdt",objective = "regression",metric = "RMSE",min_data = 1L,
                           learning_rate = 0.01,num_iterations=1000,num_leaves=3,max_depth=-1,
                           early_stopping_round=50L,cat_l2=10,skip_drop=0.5,drop_rate=0.5,
                           cat_smooth=5)
```

learning rate: The default setting is 0.1, and it is usually set between 0.05 and 0.1. Choosing a smaller learning rate can stable and better model performance.
objective: Regression presents our job is regression job.

### Others function Example
```R
result<-line_trait_mean(data=trait,trait="FTgdd",mean=env_trait,LbyE=LbyE,row=2)
MSE<-result[[1]]
ltm<-result[[2]]
mse_plotter(MSE)

Reg<-Reg(LbyE=LbyE,env_trait=env_trait)
Reg_plotter(Reg=Reg)
Mean_trait_plot(Reg,MSE)

Slope_Intercept<-Slope_Intercept(data=filtered_trait,input=env_trait, env_paras=PTT_PTR,
                  Para_Name="PTS",line="line_code",trait="PH",filter=5,
                  maxR_dap1=22, maxR_dap2=35,rounds=4)

prdM <- LOOCV(maxR_dap1=22,maxR_dap2=35,data=filtered_trait,input=env_trait,env_paras=PTT_PTR,
              Para_Name="PTS",trait="PH",line="line_code",p=1,filter=4)

prdM_plotter(prdM=prdM,data=envMeanPara,trait="PH",Para_Name="PTS");
envMeanPara_plotter(envMeanPara)
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
* GE_CV
* h2_rrBLUP
* LbyE_calculate
* LbyE_corrplot
* line_trait_mean
* LOOCV
* ltm_plotter
* Mean_trait_plot
* mse_plotter
* prdM_plotter
* Reg
* Reg_plotter
* Slope_Intercept

## Implementation notes
ggGSE is a collection of tools for cross-environmental genome-wide selection prediction that integrates most genome-wide prediction models, both parametric and non-parametric. You can input your own collected data against sample data and get the results you want directly through the built-in functions of the toolkit, which requires no additional statistical knowledge or coding skills and is somewhat user-friendly because it saves users from having to search for various tools and apply them to cross-environmental prediction.

## License
[GPL-3](https://www.gnu.org/licenses/quick-guide-gplv3.html)

## Reference
[Li X, Guo T, Wang J, et al. An integrated framework reinstating the environmental dimension for GWAS and genomic selection in crops[J]. Molecular Plant](https://www.sciencedirect.com/science/article/pii/S167420522100085X)

[Jarquín D, Crossa J, Lacaze X, et al. A reaction norm model for genomic selection using high-dimensional genomic and environmental data[J]. Theoretical and applied genetics, 2014, 127: 595-607.](https://link.springer.com/article/10.1007/s00122-013-2243-1)




