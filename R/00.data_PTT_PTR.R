#' @title PTT_PTR Dataset
#'
#' @description
#' This dataset contains information about environmental parameters related to
#' PTT (Photothermal Time) and PTR (Photothermal Ratio) and others.
#'
#' @format A data frame representing environmental parameters.
#' \describe{
#'   \item{Date}{Date of the environmental information.}
#'   \item{env_code}{Environmental code assigned to each environment.}
#'   \item{PTS}{Environmental parameter related to Photothermal spectroscopy.}
#'   \item{PTT}{Environmental parameter related to Photothermal Time.}
#'   \item{PTR}{Environmental parameter related to Photothermal Ratio.}
#'   \item{GDD}{Environmental parameter related to Growing Degree Days.}
#'   \item{Tmax}{Environmental parameter related to the maximum temperature.}
#'   \item{Tmin}{Environmental parameter related to the minimum temperature.}
#'   \item{DL}{Environmental parameter related to Day Length.}
#' }
#'
#' @source This dataset originates from the article and provides
#' environmental parameters for analysis.
#' @keywords dataset
#'
#' @usage data(PTT_PTR)
#'
#' @examples
#' data(PTT_PTR)
#'
"PTT_PTR"
