#' @title trait Dataset
#' @description
#' This dataset contains information about phenotypic traits for various
#' lines and environments.
#'
#' @format A data frame representing trait parameters.
#' \describe{
#'   \item{line_code}{Line code assigned to each line.}
#'   \item{env_code}{Environmental code assigned to each environment.}
#'   \item{env_note}{Environmental note assigned to each environment.}
#'   \item{FTdap}{Name of the phenotypic trait measured, related to FT.}
#'   \item{FTgdd}{Name of the phenotypic trait measured, related to FT.}
#'   \item{pop_code}{Other columns present additional trait-related information}
#' }
#'
#' @source This dataset originates from the article and provides
#' trait data for analysis.
#' @keywords dataset
#'
#' @usage data(trait)
#'
#' @examples
#' data(trait)
#'
"trait"
