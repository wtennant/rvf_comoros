# format_ndvi.R: Format the raw NDVI data.

# Clear the workspace.
rm(list = ls())

# Load in data manipulation libraries.
library(dplyr)
library(tidyr)

# Load in the raw data.
ndvi <- read.csv("../data/raw/ndvi.csv", header=TRUE) %>%
  mutate(DATE = as.Date(DATE)) %>%
  mutate(TIME = as.numeric(DATE - first(DATE)))

# Set a list of locations to calculate the ndvi for.
locations <- c("Anjouan",
               "Grande Comore",
               "Mayotte",
               "Moheli")

# Define a look up table for saving the information.
file_names <- c("Anjouan" = "anjouan",
                "Grande Comore" = "grande_comore",
                "Mayotte" = "mayotte",
                "Moheli" = "moheli")

# Define the Nadaraya-Watson kernel-weighted average function:
kersmooth <- function(x, y, xout=x, sigma=1)
{
  yout <- rep(0, length(xout))
  for (i in seq(1, length(xout)))
  {
    # Calculate all weights in the kernel given the current point.
    kernel_val <- exp(-(xout[i] - x)^2 / 2 / sigma / sigma)
    
    # Calculate the kernel-weighted average at point xout[i].
    yout[i] <- sum(kernel_val*y) / sum(kernel_val)
  }
  
  # Return the kernel-weighted average at given x values.
  return(data.frame("x" = xout, "y" = yout))
}

# Set up a variable to store all NDVI data.
all_ndvi <- data.frame()

# For ever location, format the NDVI data.
for (location in locations)
{
  # Calculate the kernel-weighted average NDVI for each day.
  my_data <- dplyr::filter(ndvi, ISLAND==location)
  new_data <- data.frame(TIME = seq(min(my_data$TIME), max(my_data$TIME), by=1))
  my_model <- kersmooth(my_data$TIME, my_data$NDVI, xout=new_data$TIME, sigma=21)
  new_data$NDVI <- my_model$y
  
  # Calculate the date, month and year.
  new_data <- new_data %>%
    mutate(DATE = ndvi$DATE[1] + TIME) %>%
    mutate(MONTH = as.numeric(format(DATE, format="%m"))) %>%
    mutate(YEAR = as.numeric(format(DATE, format="%Y"))) %>%
    mutate(DAY = as.numeric(format(DATE, format="%d"))) %>%
    mutate(WEEK = DAY %/% 7) %>%
    dplyr::filter(WEEK <= 3)
  
  # Average the NDVI for each month.
  new_data <- new_data %>%
    group_by(MONTH, YEAR, WEEK) %>%
    summarise(NDVI = first(NDVI),
              DATE = first(DATE)) %>%
    dplyr::select(DATE, WEEK, MONTH, YEAR, NDVI) %>%
    arrange(DATE)z
  
  # Start simulations in July 2004 and end December 2019.
  new_data <- dplyr::filter(new_data, DATE >= as.Date("2004-07-01"), DATE < as.Date("2020-07-01"))
  
  # Note the ID number of the island.
  new_data$ISLAND_ID <- which(location == locations) - 1
  new_data <- dplyr::select(new_data, ISLAND_ID, DATE, WEEK, MONTH, YEAR, NDVI)
  
  # Save the formatted version as a .csv.
  write.csv(new_data, paste0("../data/ndvi_", as.character(file_names[location]), ".csv"), row.names=FALSE)
  
  # Combine NDVI data from each island.
  all_ndvi <- bind_rows(all_ndvi, new_data)
}

# Save all the NDVI data.
write.csv(all_ndvi, "../data/ndvi.csv", row.names=FALSE)
