# prepare_post.R: Given a model fit, prepare the posterior distribution ready for scenario testing.

# Clear the workspace.
rm(list = ls())

# Load in data manipulation and visualisation libraries.
library(dplyr)
library(tidyr)
library(data.table)

# Define the burn in period.
burn_in <- 600000

# Define the number of samples to take from the posterior distribution.
n_samples <- 10000

# Define the output folder for the MCMC chains.
out_folder <- "../cpp/out/final/"

# Define the input folder for scenario testing.
in_folder <- "../cpp/in/final/"

# Get the number of chains that were run.
n_chains <- length(list.files(pattern="mcmc_[0-9]+", path=out_folder))

# Create data frame to store all the chains.
chains <- data.frame()

# Load in the MCMC chains.
for (i in seq(1, n_chains))
{
  # Load in the MCMC data.
  chain <- fread(paste0(out_folder, "mcmc_", i-1, ".csv"), header=TRUE, sep=",")
  
  # Remove the log-likelihood and prior columns from the chain data.
  chain <- select(chain, -LOG_LIKELIHOOD, -LOG_PRIOR)
  
  # Get rid of the burn in period.
  chain <- chain %>%
    filter(STEP > burn_in) %>%
    mutate(ID = (i-1)*(max(STEP) - burn_in) + (STEP - burn_in))
  
  # Store with all other chains.
  chains <- rbind(chains, chain)
}

# Take n_samples from all the chains.
sample_id <- sample(seq(min(chains$ID), max(chains$ID)), size=n_samples)
post <- chains %>%
  filter(ID %in% sample_id)

# Re-order the columns.
post <- select(post, ID, PAR_NAME, PAR_VALUE)

# Save the posterior distribution to file, ready for scenario testing.
write.csv(x=post, file=paste0(in_folder, "posterior.csv"), row.names = FALSE)

# Save a copy to the out folder as well.
write.csv(x=post, file=paste0(out_folder, "posterior.csv"), row.names = FALSE)
