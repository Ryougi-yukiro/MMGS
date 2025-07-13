#' @title Regression
#' @description
#' Calculate regression results between traits and
#' environmental means for each line.
#'
#' @param LbyE Data frame containing lines and their corresponding
#'  trait values within environments.It is from the function of LbyE_calculate.
#' @param env_trait Data frame containing environmental means for
#' different traits within environments.It is from the function
#' of enV_mean_calculate.
#' @param filter_num Minimum number of observations required to
#' perform regression. Default is 4.
#'
#' @return A data frame containing regression results for each line.
#' @export
#'
#' @examples
#' #Get Input
#' LbyE<-LbyE_calculate(trait,trait="FTgdd",env="env_code",line="line_code")
#' env_trait<-env_trait_calculate(data=trait,trait="FTgdd",env="env_code")
#' #Calculate
#' Regression<-Reg(LbyE=LbyE,env_trait=env_trait)
#' print(head(Regression))
Reg <- function(LbyE,
                env_trait,
                filter_num = NULL) {
  if (is.null(filter_num)) {
    filter_num = 4
  }
  data <- data.frame()

  order_index <-colnames(LbyE)[-1]
  env_trait_sorted <- env_trait[match(order_index, env_trait$env_code), ]
  # Loop through each line in LbyE
  for (i in 1:nrow(LbyE)) {
    # Create a data frame with meanY, Trait, and env_code
    df <- data.frame(meanY = env_trait_sorted[["mean"]],
                     Trait = as.numeric(LbyE[i, -1]),
                     env_code = colnames(LbyE)[-1])

    # Remove rows with NA trait values
    df <- df[!is.na(df$Trait),]

    # Check if the number of observations is greater than or equal to the filter_num
    if (nrow(df) >= filter_num) {
      # Add line information to the data frame
      df$line <- LbyE[i, 1]

      # Append the data frame to the overall results data frame
      data <- rbind(data, df)
    }
  }

  return(data)
}
