#' @title LbyE_Reshape
#' @description
#' Reshape trait values within specified environments and lines.
#'
#' @param data Data frame containing environmental information from env_trait_calculate.
#' @param env Column name representing different environments in the data frame.
#' @param LbyE A data frame calculated by LbyE_calculate. Column names represent Lines,
#' and row names represent Environments.
#'
#' @return Data frame containing trait values within specified environments and lines.
#' @export
#'
#' @examples
#' #Get inputs:
#' env_trait<-env_trait_calculate(data=trait,trait="FTgdd",env="env_code")
#' LbyE<-LbyE_calculate(data=trait,trait="FTgdd",env="env_code",line="line_code")
#' #Run:
#' etl<-LbyE_Reshape(data=env_trait,LbyE=LbyE,env="env_code")
#'
LbyE_Reshape <- function(data, LbyE, env) {
  env_trait <- data
  # Subset the environmental data based on the provided environments and lines
  env_trait <- env_trait[env_trait[[env]] %in% colnames(LbyE),]

  # Initialize an empty data frame to store the results
  data <- data.frame()

  # Loop through each line in LbyE
  for (i in 1:nrow(LbyE)) {
    # Create a data frame with env_code, env_order, trait, and line information
    df <- data.frame(
      env_code = colnames(LbyE)[-1],
      trait = as.numeric(LbyE[i, -1])
    )

    # Remove rows with NA trait values
    df <- df[!is.na(df$trait),]

    # Add line information to the data frame
    df$line <- LbyE[i, 1]

    # Append the data frame to the overall results data frame
    data <- rbind(data, df)
  }

  return(data)
}
