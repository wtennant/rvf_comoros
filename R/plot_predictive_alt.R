# plot_predictive_alt.R: Plot the posterior predictive distribution for the data
# given fits from scenario testing.
# These distributions are taken from simulations executed on a sample of
# the posterior distribution created by prepare_post.R

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

# Load in the serology, ndvi, simulation and chain data.
sero <- fread("../data/seroprev.csv", header=TRUE)
ndvi <- fread("../data/ndvi.csv", header=TRUE)
sim <- fread("../cpp/scenario_testing/out/scenario_free/simulation.csv", header=TRUE)

# Calculate the epidemiological year and epidemiological month in the simulated data.
sim <- sim %>%
  mutate(MONTH = ((ndvi$MONTH[1] - 1 + (TIME %/% 4)) %% 12) + 1,
         YEAR = ndvi$YEAR[1] + ((ndvi$MONTH[1] - 1 + TIME %/% 4) %/% 12),
         EPI_YEAR_START = ifelse(MONTH < 7, YEAR - 1, YEAR),
         EPI_MONTH = (MONTH + 5) %% 12 + 1)

# Get all unique date and island combinations.
unique_date_island <- dplyr::select(sero, ISLAND_ID, EPI_YEAR_START, EPI_MONTH) %>%
  unique()

# Duplicate rows for each month where data is summarised of a year.
duplicate_summaries <- filter(unique_date_island, EPI_MONTH == 13) %>%
  select(-EPI_MONTH)
duplicate_summaries <- data.frame(ISLAND_ID = rep(duplicate_summaries$ISLAND_ID, each=12),
                 EPI_YEAR_START = rep(duplicate_summaries$EPI_YEAR_START, each=12),
                 DUMMY_EPI_MONTH = rep(seq(1,12), times=nrow(duplicate_summaries)),
                 EPI_MONTH = 13)

# Add these to all unique date and island combinations.
unique_date_island <- unique_date_island %>%
  filter(EPI_MONTH != 13) %>%
  mutate(DUMMY_EPI_MONTH = EPI_MONTH) %>%
  bind_rows(duplicate_summaries) %>%
  arrange(EPI_YEAR_START, EPI_MONTH, DUMMY_EPI_MONTH, ISLAND_ID)

# Filter out all relevant simulation data.
relevant_sim <- inner_join(unique_date_island, sim, by=c("ISLAND_ID", "EPI_YEAR_START", "DUMMY_EPI_MONTH"="EPI_MONTH"))

# Calculate summaries for each age, across all ages, as well as adults.
each_age_sum <- group_by(relevant_sim, SAMPLE_ID, ISLAND_ID, EPI_YEAR_START, EPI_MONTH, AGE) %>%
  summarise(SERO_PREV = sum(R) / sum(S + E + I + R))
all_age_sum <- group_by(relevant_sim, SAMPLE_ID, ISLAND_ID, EPI_YEAR_START, EPI_MONTH) %>%
  summarise(SERO_PREV = sum(R) / sum(S + E + I + R)) %>%
  mutate(AGE = -1)
adult_sum <- filter(relevant_sim, AGE != 0) %>%
  group_by(SAMPLE_ID, ISLAND_ID, EPI_YEAR_START, EPI_MONTH,) %>%
  summarise(SERO_PREV = sum(R) / sum(S + E + I + R)) %>%
  mutate(AGE = -2)
young_sum <- filter(each_age_sum, AGE == 0) %>%
  mutate(AGE = -3)

# Union all simulation summaries.
sim_sero <- bind_rows(each_age_sum, all_age_sum) %>%
  bind_rows(adult_sum) %>%
  bind_rows(young_sum)

# Label the ID of each data entry and calculate seroprevalence.
sero <- sero %>%
  mutate(SERO_PREV = N_POSITIVE / N_TESTED) %>%
  mutate(LOWER_CI = binom.confint(N_POSITIVE, N_TESTED, method="exact")$lower) %>%
  mutate(UPPER_CI = binom.confint(N_POSITIVE, N_TESTED, method="exact")$upper)

# Calculate the year and month of empirical and simulated data.
sero <- sero %>%
  mutate(MONTH = ifelse(EPI_MONTH != 13, (EPI_MONTH - 7) %% 12 + 1, 13),
         YEAR = ifelse(MONTH >= 7, EPI_YEAR_START, EPI_YEAR_START + 1))
sim_sero <- sim_sero %>%
  mutate(MONTH = ifelse(EPI_MONTH != 13, (EPI_MONTH - 7) %% 12 + 1, 13),
         YEAR = ifelse(MONTH >= 7, EPI_YEAR_START, EPI_YEAR_START + 1))

# Rename AGE category of simulated serological data.
sim_sero <- sim_sero %>%
  rename(AGE_GROUP = AGE)

# Define the age labels.
age_labels <- paste(seq(0, 9), seq(1, 10), sep="–")
age_labels[10] <- "9+"

# Define the labels for each facet.
# Linking ISLAND_IDs to island names.
island_name <- c("1" = "Grande Comore",
                 "3" = "Mohéli",
                 "0" = "Anjouan",
                 "2" = "Mayotte")

# Strings for calendar years and months.
year_name <- as.character(seq(2004, 2020))
names(year_name) <- year_name
month_name <- c(month.abb[seq(1,12)], "All")
names(month_name) <- as.character(seq(1,13))

# Create all combinations of months, year and island.
facet_df <- expand.grid(MONTH = names(month_name),
                        YEAR = as.integer(names(year_name)),
                        ISLAND_ID = names(island_name))

# Create the labels associated with each combintion.
facet_df <- facet_df %>%
  mutate(ISLAND_NAME = island_name[ISLAND_ID],
         YEAR_NAME = year_name[as.character(YEAR)],
         MONTH_NAME = month_name[MONTH]) %>%
  mutate(LABEL = paste(ISLAND_ID, YEAR, MONTH, sep="-"),
         LABEL_NAME = paste0(MONTH_NAME, " ", YEAR)) %>%
  mutate(LABEL_NAME = ifelse(MONTH==13, paste0(YEAR, "—", YEAR + 1), LABEL_NAME))

# Extract all label information.
facet_label <- facet_df$LABEL_NAME
names(facet_label) <- facet_df$LABEL

# Change the group ID.
sero <- sero %>%
  mutate(GROUP_ID = as.numeric(factor(paste(AGE_GROUP, ISLAND_ID, YEAR, MONTH)))) %>%
  mutate(LABEL = factor(paste(ISLAND_ID, YEAR, MONTH, sep="-"), levels=names(facet_label))) %>%
  mutate(LABEL = factor(LABEL, levels=unique(LABEL)))
sim_sero <- sim_sero %>%
  mutate(GROUP_ID = as.numeric(factor(paste(AGE_GROUP, ISLAND_ID, YEAR, MONTH)))) %>%
  mutate(LABEL = factor(paste(ISLAND_ID, YEAR, MONTH, sep="-"), levels=names(facet_label))) %>%
  mutate(LABEL = factor(LABEL, levels=unique(LABEL)))

# Re-order age groups.
sero <- sero %>%
  mutate(AGE_GROUP = ifelse(AGE_GROUP == -1, -4, AGE_GROUP)) %>%
  mutate(AGE_GROUP = ifelse((AGE_GROUP == -2) | (AGE_GROUP == -3), AGE_GROUP + 1, AGE_GROUP))
sim_sero <- sim_sero %>%
  mutate(AGE_GROUP = ifelse(AGE_GROUP == -1, -4, AGE_GROUP)) %>%
  mutate(AGE_GROUP = ifelse((AGE_GROUP == -2) | (AGE_GROUP == -3), AGE_GROUP + 1, AGE_GROUP)) %>%
  mutate(AGE_GROUP = ifelse(AGE_GROUP >= 0, AGE_GROUP + 1, AGE_GROUP))

# Visualise the predictive posterior distributions.
gg_pred <- list()
for (island_id in seq(0,3))
{
  # Filter the serological data for the island.
  island_sim_sero <- filter(sim_sero, ISLAND_ID == island_id)
  island_sero <- filter(sero, ISLAND_ID == island_id)
  
  # If Anjouan, add one row below,
  if (island_id == 0) {
    island_sim_sero$LABEL <- factor(island_sim_sero$LABEL, levels=c(as.character(unique(island_sim_sero$LABEL)), "4", "5", "6"))
    island_sero$LABEL <- factor(island_sero$LABEL, levels=c(as.character(unique(island_sero$LABEL)), "4", "5", "6"))
  } else if (island_id == 3) {  # If Mohéli, add two rows above.
    island_sim_sero$LABEL <- factor(island_sim_sero$LABEL, levels=c("4", "5", "6", "7", "8", "9", as.character(unique(island_sim_sero$LABEL))))
    island_sero$LABEL <- factor(island_sero$LABEL, levels=c("4", "5", "6", "7", "8", "9", as.character(unique(island_sero$LABEL))))
  } else if (island_id == 2) {  # For Mayotte, add one row above.
    island_sim_sero$LABEL <- factor(island_sim_sero$LABEL, levels=c("4", "5", "6", as.character(unique(island_sim_sero$LABEL))))
    island_sero$LABEL <- factor(island_sero$LABEL, levels=c("4", "5", "6", as.character(unique(island_sero$LABEL))))
  } else {   # For Grande Comore, don't add any more facets!
    island_sim_sero$LABEL <- factor(island_sim_sero$LABEL, levels=c(as.character(unique(island_sim_sero$LABEL))))
    island_sero$LABEL <- factor(island_sero$LABEL, levels=c(as.character(unique(island_sero$LABEL))))
  }
  
  # Make the plot!
  gg_pred[[island_id + 1]] <- ggplot(island_sim_sero, aes(x=AGE_GROUP, y=SERO_PREV, group=GROUP_ID)) +
    geom_errorbar(data=island_sero, aes(ymin=LOWER_CI, ymax=UPPER_CI, alpha="DATA"), width=1, colour="grey75") +
    geom_point(data=island_sero, aes(alpha="DATA"), colour="grey60") +
    geom_violin(scale="width", aes(fill=as.factor(ISLAND_ID), colour=as.factor(ISLAND_ID)), alpha=0.75) +
    geom_vline(xintercept=0, linetype="dashed", colour="grey50") +
    geom_vline(xintercept=-3, linetype="dashed", colour="grey50") +
    scale_x_continuous(breaks=c(-4, -2, -1, seq(1,10)),
                       labels=c("Unknown", "Young", "Adult", age_labels)) +
    scale_y_continuous(breaks=c(0, 0.5, 1.0)) +
    scale_colour_brewer("Simulated data", palette="Set1",
                        limits=names(island_name),
                        labels=island_name,
                        guide=FALSE) +
    scale_fill_brewer("Simulated data", palette="Pastel1",
                      labels=island_name,
                      limits=names(island_name),
                      guide=FALSE) +
    scale_alpha_manual("Empirical evidence", breaks="DATA", values=1, label="Mean with 95%\nconfidence interval", guide=FALSE) +
    xlab("Age group") +
    ylab("IgG seroprevalence") +
    facet_wrap(~LABEL, labeller=labeller(LABEL=facet_label), ncol=3, drop=FALSE) +
    gg_theme +
    theme(panel.border=element_rect(colour="black", fill=NA)) +
    theme(axis.text.x = element_text(angle=90, size=10, hjust=1, vjust=0.5)) +
    theme(legend.position="right") +
    theme(strip.text = element_text(size=12)) +
    #guides(colour=guide_legend(title.hjust=0, title.position="top", order=0),
    #       fill=guide_legend(order=0), alpha=guide_legend(order=1)) +
    ggtitle(as.character(island_name[as.character(island_id)]))
}

## Function to extract legend
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

# Construct the legend
gg_legend <- ggplot(island_sim_sero, aes(x=AGE_GROUP, y=SERO_PREV, group=GROUP_ID)) +
  geom_errorbar(data=island_sero, aes(ymin=LOWER_CI, ymax=UPPER_CI, alpha="DATA"), width=1, colour="grey75") +
  geom_point(data=island_sero, aes(alpha="DATA"), colour="grey60") +
  geom_violin(scale="width", aes(fill=as.factor(ISLAND_ID), colour=as.factor(ISLAND_ID)), alpha=0.75) +
  scale_colour_brewer("Simulated data", palette="Set1",
                      limits=names(island_name),
                      labels=island_name) +
  scale_fill_brewer("Simulated data", palette="Pastel1",
                    labels=island_name,
                    limits=names(island_name)) +
  scale_alpha_manual("Empirical evidence", breaks="DATA", values=1, label="Mean with 95%\nconfidence interval") +
  facet_wrap(~LABEL, labeller=labeller(LABEL=facet_label), ncol=3, drop=FALSE) +
  gg_theme +
  guides(colour=guide_legend(title.hjust=0, title.position="top", order=0),
         fill=guide_legend(order=0), alpha=guide_legend(order=1, title.position="top"))

# Convert to gg plot object.
gg_legend <- ggpubr::as_ggplot(g_legend(gg_legend))

# Arrange all in a grid.
p1 <- egg::ggarrange(gg_pred[[2]], gg_pred[[1]], ncol=2, draw=FALSE)
p2 <- egg::ggarrange(gg_pred[[4]], gg_pred[[3]], ncol=2, draw=FALSE)
p3 <- ggpubr::ggarrange(gg_legend, p1, p2, nrow=3, heights=c(1, 5, 5))

# Save the figure and crop.
file_name <- "model_predictive"
ggsave(filename=paste0(file_name, ".svg"), plot=p3, width=16*1.05, height=16*1.2, dpi=600, scale=1.0,
       path="figures")
ggsave(filename=paste0(file_name, ".pdf"), plot=p3, width=16*1.05, height=16*1.2, dpi=600, scale=1.0,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, ".pdf figures/", file_name, ".pdf\""))
