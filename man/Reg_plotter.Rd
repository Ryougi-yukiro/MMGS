% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03.Regression_plotter.R
\name{Reg_plotter}
\alias{Reg_plotter}
\title{Reg_plotter}
\usage{
Reg_plotter(Reg = Reg, size = NULL, method = NULL, color = NULL, alpha = NULL)
}
\arguments{
\item{Reg}{Input data frame containing regression results.
param x Variable to be plotted on the x-axis.
param y Variable to be plotted on the y-axis.
param group Variable used for grouping the data.}

\item{size}{Size of points in the plot.}

\item{method}{Smoothing method for the regression line.}

\item{color}{Vector containing colors for regression line and points.}

\item{alpha}{Vector containing alpha values for regression line and points.}
}
\value{
A plot for the regression relationship between mean and trait.
}
\description{
Generate a regression plot between the mean and trait,
showing the relationship between these two variables.
}
\examples{
#' #Get Input
LbyE<-LbyE_calculate(trait,trait="FTgdd",env="env_code",line ="line_code")
env_trait<-env_trait_calculate(data=trait,trait ="FTgdd",env ="env_code")
Regression<-Reg(LbyE = LbyE, env_trait = env_trait)
Reg_plotter(Reg = Regression)

}
