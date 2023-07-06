# ggGSE
A  devlopment R packages used for predicting traits from different envs


## Description

A brief description of your R package, highlighting its purpose and main functionalities.

## Installation

You can install the package from CRAN using the following command:
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

#Get Envs mean within env paras
envMeanPara<-envMeanPara(data=env_trait, env_paras=PTT_PTR, maxR_dap1=18,maxR_dap2=43, Paras=Paras)
```

### Step.3 CV
Users can customize the model they need, the function uses the norm reaction by default, the given environment parameters can be obtained from the previous results , fold number represents the number of folds, reshuffle represents the number of repetitions. In RM.G mode, the available models are rrBLUP, LASSO,EN,RR,BA,BB,BC,BL,BRR,RKHS,MKRKHS,SVM,RF and LightGBM.
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

#Others function
```R
result<-line_trait_mean(data=trait,trait="FTgdd",mean=env_trait,LbyE=LbyE,row=2)
MSE<-result[[1]]
ltm<-result[[2]]
mse_plotter(MSE)
```

## Vignettes
If your package includes vignettes, mention them here and provide instructions on how to access and utilize them.

## Issues
If you encounter any issues or have suggestions for improvements, please open an issue on the GitHub repository.


## License
Include information about the license under which your package is released.


## Citation
Li X, Guo T, Wang J, et al. An integrated framework reinstating the environmental dimension for GWAS and genomic selection in crops[J]. Molecular Plant<https://www.sciencedirect.com/science/article/pii/S167420522100085X>

Jarquín D, Crossa J, Lacaze X, et al. A reaction norm model for genomic selection using high-dimensional genomic and environmental data[J]. Theoretical and applied genetics, 2014, 127: 595-607.
![image](https://github.com/Ryougi-yukiro/ggGSE/assets/64737958/2c716e1d-d980-4037-bcf3-21af9b5978c6)




