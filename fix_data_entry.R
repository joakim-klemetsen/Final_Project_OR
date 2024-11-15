# load required packages
library(tidyverse)

# load data
data <- read.csv("data/ProjectData.csv", sep = ";", stringsAsFactors = FALSE)

# function fixing issues with thousand separator in data entry
number_fix <- function(variable) {
  result <- numeric(length(variable)) 
  for (i in seq_along(variable)) {
    split <- str_split(variable[i], "\\.", simplify = TRUE) %>% as.data.frame()
    if (ncol(split) > 2) {
      result[i] <- 1e3 * as.numeric(split[1]) + as.numeric(split[2]) + (1 / 1e3) * as.numeric(split[3])
    } else if (ncol(split) == 2) {
      result[i] <- as.numeric(split[1]) + (1 / 1e3) * as.numeric(split[2])
    } else {
      result[i] <- as.numeric(split[1])
    }
  }
  return(result)
}

# apply the function
data$Price <- sapply(data$Price, function(x) number_fix(x))

# write the modified data frame to a new CSV file
write.csv(data, "data/ModifiedProjectData.csv", row.names = FALSE)
