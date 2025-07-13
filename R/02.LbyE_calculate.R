#' @title LbyE_calculate
#' @description
#' Calculate the lines in all environments.
#'
#' @param env Column name representing different environments
#' in the input data frame.
#' @param line Column name representing different lines in the input data frame.
#' @param data Input data frame containing trait values, line information,
#' and environmental information.
#' @param trait The trait for which the lines within different environments
#' will be calculated.
#' @return A data frame containing line values within different environments
#' for the specified trait.
#' @export
#'
#' @examples
#' # Calculate lines within different environments for the trait "FTgdd"
#' LbyE_calculate(trait,trait="FTgdd",env="env_code",line="line_code")

LbyE_calculate<-function(data,trait,env,line){
  # Get unique environment and line codes
  envs <- unique(as.vector(data[[env]]))
  lines <- unique(as.vector(data[[line]]))
  # Create a data frame to store the results
  LbyE <- data.frame(line_code = lines)

  # Find the column indices for trait and line in the data frame
  t <- match(trait, colnames(data))
  t2 <- match(line, colnames(data))
  # Loop through each environment
  for (e_l in 1:length(envs)) {
    e <- envs[e_l]

    # Subset data for the current environment
    e_l_trait <- subset(data, data[[env]] == e)

    # Rename the trait column to the current environment code
    nonNAs <- length(which(!is.na(e_l_trait[, t])))
    colnames(e_l_trait)[t] <- e

    # Merge the subsetted data with the results data frame
    LbyE <- merge(LbyE, e_l_trait[, c(t2, t)], all.x = TRUE)
  }

  return(LbyE)
}
