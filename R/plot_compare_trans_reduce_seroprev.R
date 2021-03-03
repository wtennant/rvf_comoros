# plot_compare_trans_reduce_seroprev.R: Plot the effects of different transmission
# reduction scenarios on seroprevalence over time.

# Clear the workspace.
rm(list = ls())

# Load in the data visualisation and manipulation libraries.
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)

# Load in our favourite ggplot theme.
source("gg_theme.R")

# Define the directory with the scenarios in.
folder <- "../cpp/transmission_reduce/out/"

# Define a look up table with the scenario name and the
# sub-folder its contained in
scenarios <- c("Full transmission" = "no_control/",
               "Grande Comore (10%)" = "gc_10/",
               "Grande Comore (20%)" = "gc_20/",
               "Grande Comore (30%)" = "gc_30/",
               "Mohéli (10%)" = "moh_10/",
               "Mohéli (20%)" = "moh_20/",
               "Mohéli (30%)" = "moh_30/",
               "Anjouan (10%)" = "anj_10/",
               "Anjouan (20%)" = "anj_20/",
               "Anjouan (30%)" = "anj_30/",
               "Mayotte (10%)" = "may_10/",
               "Mayotte (20%)" = "may_20/",
               "Mayotte (30%)" = "may_30/")

# Load in the NDVI data.
ndvi <- fread("../data/ndvi.csv", header=TRUE, sep=",")

# Use to get dates corresponding to the times in the simulation.
dates <- ndvi %>%
  mutate(DATE = as.Date(DATE)) %>%
  filter(ISLAND_ID == 0)

# The last date we need is one week after the NDVI data finishes.
final_date <- data.frame("ISLAND_ID" = 0,
                   "WEEK" = (dates$WEEK[nrow(dates)] + 1) %% 4,
                   "MONTH" = (dates$MONTH[nrow(dates)] + (dates$WEEK[nrow(dates)] + 1) %/% 4 - 1) %% 12 + 1,
                   "YEAR" = dates$YEAR[nrow(dates)] + ((dates$MONTH[nrow(dates)] + (dates$WEEK[nrow(dates)] + 1) %/% 4 - 1) %/% 12)) %>%
  mutate(DAY = WEEK*7 + 1) %>%
  mutate(DATE = as.Date(paste(YEAR, sprintf("%.2d", MONTH), sprintf("%.2d", DAY)), format="%Y %m %d"))

# Append the final date and just extract dates and simulation times.
dates <- dates %>%
  bind_rows(final_date) %>%
  mutate(TIME = 0:(n()-1)) %>%
  ungroup() %>%
  select(TIME, DATE)

# Create space to store the seroprev of all scenarios.
sum_scenario <- data.frame()

# For each scenario, calculate the total number of infections
# on each island for each simulation.
for (i in seq(1, length(scenarios)))
{
  # Load in the simulation data.
  sim <- fread(paste0(folder, scenarios[i], "simulation.csv"), header=TRUE)
  
  # Calculate the seroprevalence at each time point.
  sim <- group_by(sim, SAMPLE_ID, ISLAND_ID, TIME) %>%
    summarise(SERO_PREV = sum(R) / sum(S + E + I + R))
  
  # Calculate the median, 50% credible interval and 95% credible interval
  # of serprevalence at each time poitn for each island.
  sim <- sim %>%
    group_by(ISLAND_ID, TIME) %>%
    summarise(Q_500 = quantile(SERO_PREV, probs=0.5),
              Q_250 = quantile(SERO_PREV, probs=0.25),
              Q_750 = quantile(SERO_PREV, probs=0.75),
              Q_975 = quantile(SERO_PREV, probs=0.975),
              Q_025 = quantile(SERO_PREV, probs=0.025))
  
  # Save the scenario ID and store in the overall data frame
  # of summary metrics.
  sum_scenario <- sim %>%
    mutate(SCENARIO_ID = i) %>%
    bind_rows(sum_scenario, .)
}

# Join with the corresponding dates.
sum_scenario <- left_join(sum_scenario, dates, by=c("TIME"))

# Define the names of each island based on their ID number.
island_labels <- c("0" = "Anjouan",
                  "1" = "Grande Comore",
                  "2" = "Mayotte",
                  "3" = "Mohéli")

# Define the facet labels.
facet_labels <- names(scenarios)
names(facet_labels) <- paste(seq(1, length(scenarios)))

# Change ISLAND_ID and SCENARIO_ID to a factor.
sum_scenario$SCENARIO_ID <- as.factor(sum_scenario$SCENARIO_ID)
sum_scenario$ISLAND_ID <- factor(sum_scenario$ISLAND_ID, levels=c(1, 3, 0, 2))

# Plot seroprevalence over time.
gg_sp <- ggplot(sum_scenario, aes(x=DATE)) +
  geom_ribbon(aes(ymin=Q_025, ymax=Q_975, fill=ISLAND_ID, alpha="95_CI")) +
  geom_ribbon(aes(ymin=Q_250, ymax=Q_750, fill=ISLAND_ID, alpha="50_CI")) +
  geom_line(aes(y=Q_500, colour=ISLAND_ID, alpha="MED"), size=1) +
  geom_vline(xintercept=as.Date("2015-07-01"), linetype="dashed") +
  scale_x_date(date_breaks="1 year", date_labels="%Y") +
  scale_y_continuous(limits=c(0,1), labels=scales::percent) +
  scale_colour_brewer("Island", palette="Set1",
                      labels=island_labels[levels(sum_scenario$ISLAND_ID)]) +
  scale_fill_brewer("Island", palette="Set1",
                    labels=island_labels) +
  scale_alpha_manual("Posterior distribution",
                     breaks=c("MED", "50_CI", "95_CI"),
                     labels=c("MED" = "Median", "50_CI" = "50% credible interval", "95_CI" = "95% credible interval"),
                     values=c("MED" = 0.75, "50_CI" = 0.3, "95_CI" = 0.2)) +
  xlab("Calendar year") +
  ylab("Percentage of animals seropositive") +
  facet_wrap(~SCENARIO_ID, ncol=1, labeller=labeller(SCENARIO_ID=facet_labels)) +
  gg_theme +
  theme(panel.border=element_rect(colour="black", fill=NA)) +
  guides(colour=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                            override.aes=list(alpha=0.5)),
        alpha=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                           override.aes=list(linetype=c(1,0,0),
                                             alpha=c(1, 0.4, 0.2),
                                             fill=c(NA, "grey25", "grey25"))))
print(gg_sp)

# Save figure as an svg and pdf.
file_name <- "scenario_vec_control_seroprev"
ggsave(filename=paste0(file_name, ".svg"), plot=gg_sp, width=16, height=2*16, dpi=600, scale=0.9,
       path="figures")
ggsave(filename=paste0(file_name, ".pdf"), plot=gg_sp, width=16, height=2*16, dpi=600, scale=0.9,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, ".pdf figures/", file_name, ".pdf\""))
