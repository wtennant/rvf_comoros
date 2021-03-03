# plot_predictive.R: Plot the posterior predictive distribution for the data.
# This distribution is calculated directly from a sample of simulations
# saved during the MCMC algorithm.

# Clear the workspace.
rm(list = ls())

# Load data manipulation and visualisation libraries.
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(binom)

# Load our favourite ggplot theme.
source("gg_theme.R")

# Set the burn in period for the chains.
burn_in <- 500000

# Load in the serology, ndvi, simulation and chain data.
sero <- fread("../data/seroprev.csv", header=TRUE)
ndvi <- fread("../data/ndvi.csv", header=TRUE)
sim <- fread("../cpp/out/simulation_1.csv", header=TRUE)

# Calculate the epidemiological year and epidemiological month in the simulated data.
sim <- sim %>%
  filter(ITER >= burn_in) %>%
  mutate(WEEK = TIME %% 4,
         MONTH = ((ndvi$MONTH[1] - 1 + (TIME %/% 4)) %% 12) + 1,
         YEAR = ndvi$YEAR[1] + ((ndvi$MONTH[1] - 1 + TIME %/% 4) %/% 12),
         EPI_YEAR_START = ifelse(MONTH < 7, YEAR - 1, YEAR),
         EPI_MONTH = (MONTH + 5) %% 12 + 1)

# Create a space to store the seroprevalence of the simulations at each data point.
sim_sero <- data.frame()

# Get all unique date and island combinations.
unique_date_island <- dplyr::select(sero, ISLAND_ID, EPI_YEAR_START, EPI_YEAR_END, EPI_MONTH) %>%
  unique()

# Expand these combinations across all possible age categories.
dummy_sero <- unique_date_island[rep(1:nrow(unique_date_island), each=13),]
dummy_sero$AGE_GROUP <- rep(c(-3, -2, -1, seq(1, 10)), times=nrow(unique_date_island))

# Set up a text progress bar.
cat("\nCalculating seroprevalence...\n")
pb <- txtProgressBar(max=nrow(dummy_sero), style=3)

# For each date, island and age group, calculate
# the seroprevalence from the simulations.
for (i in seq(1, nrow(dummy_sero)))
{
  #...calculate the ages to calculate seroprevalence for.
  if (dummy_sero$AGE_GROUP[i] == -1) {
    ages <- seq(0, 9)
  } else if (dummy_sero$AGE_GROUP[i] == -2) {
    ages <- seq(1, 9)
  } else if (dummy_sero$AGE_GROUP[i] == -3) {
    ages <- 0
  } else {
    ages <- dummy_sero$AGE_GROUP[i] - 1
  }
  
  # Calculate the epidemiological months to calculate seroprevalence for.
  if (dummy_sero$EPI_MONTH[i] == 13) {
    months <- seq(1, 12)
  } else {
    months <- dummy_sero$EPI_MONTH[i]
  }
  
  # Filter out the simulation data to correspond to the period of the empirical data.
  one_sim <- filter(sim,
                     ISLAND_ID == dummy_sero$ISLAND_ID[i],
                     EPI_YEAR_START == dummy_sero$EPI_YEAR_START[i],
                     AGE %in% ages,
                     EPI_MONTH %in% months)
  
  # Calculate the mean seroprevalence.
  sero_prev <- one_sim %>%
    group_by(ITER) %>%
    summarise(SERO_PREV = sum(R) / sum(S + E + I + R)) %>%
    mutate(DATA_ID = i,
           EPI_YEAR_START = dummy_sero$EPI_YEAR_START[i],
           EPI_MONTH = dummy_sero$EPI_MONTH[i],
           ISLAND_ID = dummy_sero$ISLAND_ID[i],
           AGE_GROUP = dummy_sero$AGE_GROUP[i])
  
  # Append to the data frame to the seroprevalence at all data points.
  sim_sero <- bind_rows(sim_sero, sero_prev)
  
  # Update progress bar.
  setTxtProgressBar(pb, i)
}

# Close the text progress bar.
close(pb)

# Label the ID of each data entry and calculate seroprevalence.
sero <- sero %>%
  mutate(DATA_ID = 1:n()) %>%
  mutate(SERO_PREV = N_POSITIVE / N_TESTED) %>%
  mutate(LOWER_CI = binom.confint(N_POSITIVE, N_TESTED, method="exact")$lower) %>%
  mutate(UPPER_CI = binom.confint(N_POSITIVE, N_TESTED, method="exact")$upper)

# Define the age labels.
age_labels <- paste(seq(0, 9), seq(1, 10), sep="–")
age_labels[10] <- "9+"

# Define the labels for each facet.
island_name <- c("1" = "Grande Comore",
                 "3" = "Mohéli",
                 "0" = "Anjouan",
                 "2" = "Mayotte")
epi_year_name <- as.character(seq(2004, 2019))
epi_year_name <- paste(epi_year_name, substr(as.character(seq(2005, 2020)), 3, 4), sep="–")
names(epi_year_name) <- as.character(seq(2004, 2019))
month_name <- c(month.abb[c(seq(7,12), seq(1,6))], "All")
names(month_name) <- as.character(seq(1,13))
facet_df <- expand.grid(EPI_MONTH = names(month_name),
                        EPI_YEAR_START = names(epi_year_name),
                        ISLAND_ID = names(island_name))
facet_df <- facet_df %>%
  mutate(ISLAND_NAME = island_name[ISLAND_ID],
         EPI_YEAR_START_NAME = epi_year_name[EPI_YEAR_START],
         EPI_MONTH_NAME = month_name[EPI_MONTH]) %>%
  mutate(LABEL = paste(ISLAND_ID, EPI_YEAR_START, EPI_MONTH, sep="-"),
         LABEL_NAME = paste(ISLAND_NAME, EPI_YEAR_START_NAME, EPI_MONTH_NAME, sep=", "))
facet_label <- facet_df$LABEL_NAME
names(facet_label) <- facet_df$LABEL

# Change the group ID.
sero <- sero %>%
  mutate(GROUP_ID = as.numeric(factor(paste(AGE_GROUP, ISLAND_ID, EPI_YEAR_START, EPI_MONTH)))) %>%
  mutate(LABEL = factor(paste(ISLAND_ID, EPI_YEAR_START, EPI_MONTH, sep="-"), levels=names(facet_label)))
sim_sero <- sim_sero %>%
  mutate(GROUP_ID = as.numeric(factor(paste(AGE_GROUP, ISLAND_ID, EPI_YEAR_START, EPI_MONTH)))) %>%
  mutate(LABEL = factor(paste(ISLAND_ID, EPI_YEAR_START, EPI_MONTH, sep="-"), levels=names(facet_label)))

# Re-order age groups.
sero <- sero %>%
  mutate(AGE_GROUP = ifelse(AGE_GROUP == -1, -4, AGE_GROUP)) %>%
  mutate(AGE_GROUP = ifelse((AGE_GROUP == -2) | (AGE_GROUP == -3), AGE_GROUP + 1, AGE_GROUP))
sim_sero <- sim_sero %>%
  mutate(AGE_GROUP = ifelse(AGE_GROUP == -1, -4, AGE_GROUP)) %>%
  mutate(AGE_GROUP = ifelse((AGE_GROUP == -2) | (AGE_GROUP == -3), AGE_GROUP + 1, AGE_GROUP))

# Visualise the predictive posterior distributions.
gg_pred <- ggplot(sim_sero, aes(x=AGE_GROUP, y=SERO_PREV, group=GROUP_ID)) +
  geom_errorbar(data=sero, aes(ymin=LOWER_CI, ymax=UPPER_CI, alpha="DATA"), width=1, colour="grey75") +
  geom_point(data=sero, aes(alpha="DATA"), colour="grey60") +
  geom_violin(scale="width", aes(fill=as.factor(ISLAND_ID), colour=as.factor(ISLAND_ID))) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey50") +
  geom_vline(xintercept=-3, linetype="dashed", colour="grey50") +
  scale_x_continuous(breaks=c(-4, -2, -1, seq(1,10)),
                     labels=c("Unknown", "Young", "Adult", age_labels)) +
  scale_y_continuous(breaks=c(0, 0.5, 1.0)) +
  scale_colour_brewer("Simulated", palette="Set1",
                      limits=names(island_name),
                       labels=island_name) +
  scale_fill_brewer("Simulated", palette="Pastel1",
                    labels=island_name,
                    limits=names(island_name)) +
  scale_alpha_manual("Data", breaks="DATA", values=1, label="Mean with 95%\nconfidence interval") +
  xlab("Age group") +
  ylab("Seroprevalence") +
  facet_wrap(~LABEL, labeller=labeller(LABEL=facet_label)) +
  gg_theme +
  theme(panel.border=element_rect(colour="black", fill=NA)) +
  theme(axis.text.x = element_text(angle=60, size=8)) +
  theme(legend.position="right") +
  guides(colour=guide_legend(title.hjust=0, title.position="top", order=0),
         fill=guide_legend(order=0), alpha=guide_legend(order=1))
print(gg_pred)

# Save the figure and crop.
file_name <- "model_predictive"
ggsave(filename=paste0(file_name, ".svg"), plot=gg_pred, width=16, height=9, dpi=600, scale=1.3,
       path="figures")
ggsave(filename=paste0(file_name, ".pdf"), plot=gg_pred, width=16, height=9, dpi=600, scale=1.3,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, ".pdf figures/", file_name, ".pdf\""))
