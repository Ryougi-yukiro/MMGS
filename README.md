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

#Step.1 Load packages and data
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
# Step.2 estimate h2 within different Envs
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



## Vignettes
If your package includes vignettes, mention them here and provide instructions on how to access and utilize them.

## Issues
If you encounter any issues or have suggestions for improvements, please open an issue on the GitHub repository.


## License
Include information about the license under which your package is released.


## Citation
If you use this package in your research or work, please consider citing it. Provide the relevant citation details here.
