#' @title line_trait_mean
#' @description
#' Calculate the mean of residuals for each line after fitting a linear model.
#'
#' @param data Input data frame containing trait values, line information,
#' and environmental information.
#' @param LbyE Data frame calculated by LbyE_calculate,
#' containing line values within different environments.
#' @param mean Data frame containing mean trait values
#' within different environments.
#' This data frame get from the function of env_mean_calculate.
#' @param trait The trait for which the residuals will be calculated.
#' @param row The minimum number of rows required for
#' fitting a linear model for each line.
#'
#' @return A list containing two elements:
#'   1. data: Data frame with columns env_code, errors (mean of residuals),
#'    and additional columns from env_trait.
#'   2. lm_residuals: Data frame containing the residuals squared for
#'   each line and environment.
#' @importFrom dplyr left_join
#' @export
#'
#' @examples
#' LbyE<-LbyE_calculate(trait,trait="FTgdd",env="env_code",line="line_code")
#' env_trait<-env_trait_calculate(data=trait,trait="FTgdd",env="env_code")
#' #Run
#' result<-line_trait_mean(data=trait,LbyE=LbyE,
#'                        mean = env_trait,trait="FTgdd",row=2)
#' print(result[[1]])
#'
#'
line_trait_mean <- function(data, LbyE, mean, trait, row) {
  # Get unique line codes from LbyE
  line_codes <- as.vector(unique(LbyE[, 1]))

  # Initialize a data frame to store residuals
  lm_residuals <- data.frame(env_code = mean[, 1])

  # Loop through each line
  for (i in line_codes) {
    # Subset data for the current line
    line_data_0 <- subset(data, data[["line_code"]] == i)

    # Check if there are enough rows for fitting a linear model
    if (nrow(line_data_0) > row) {
      # Merge data with env_trait and line_data_0
      line_data_0 <- merge(mean, line_data_0)

      # Fit a linear model and calculate residuals squared
      line_lm <- stats::lm(line_data_0[[trait]] ~ line_data_0$mean)
      df1 <- data.frame(env_code = line_data_0$env_code,
                      line_code = round((line_lm$residuals)^2, 3))

      # Rename the column with the current line code
      colnames(df1)[2] <- i

      # Merge residuals data with lm_residuals
      lm_residuals <- merge(lm_residuals, df1, all.x = TRUE)
    }
  }

  # Create a data frame with env_code, errors (mean of residuals), and additional columns from env_trait
  data <- data.frame(env_code = lm_residuals[, 1],
                     errors = rowMeans(lm_residuals[, -1], na.rm = TRUE)) %>%
    left_join(mean)

  return(list(MSE=data, lmlm_residuals=lm_residuals))
}

