# plot_compare_scenario_seroprev.R: Compares the seroprevalence of 
# all isolation scenarios against the full trade network.

# Clear the workspace.
rm(list = ls())

# Load in the data visualisation and manipulation libraries.
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)

# The date window to display the time series over.
first_date <- as.Date("2004-07-01")
last_date <- as.Date("2015-07-01")

# Load in our favourite ggplot theme.
source("gg_theme.R")

# Define the directory with the scenarios in.
folder <- "../cpp/scenario_testing/out/"

# Define a look up table with the scenario name and the
# sub-folder its contained in
scenarios <- c("Full trade network" = "scenario_free/",
               "Grande Comore imports and exports reduced by 100%" = "isolate_grande_comore/",
               "Mohéli imports and exports reduced by 100%" = "isolate_moheli/",
               "Anjouan imports and exports reduced by 100%" = "isolate_anjouan/",
               "Mayotte imports and exports reduced by 100%" = "isolate_mayotte/")

# Load in the NDVI data.
ndvi <- fread("../data/ndvi.csv", header=TRUE, sep=",")

# Use to get dates corresponding to the times in the simulation.
dates <- fread("../data/ndvi.csv", header=TRUE, sep=",") %>%
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
all_sims <- data.frame()

# For each scenario, calculate the total number of infections
# on each island for each simulation.
for (i in seq(1, length(scenarios)))
{
  # Load in the simulation data.
  sim <- fread(paste0(folder, scenarios[i], "simulation.csv"), header=TRUE, data.table=FALSE)
  
  # Calculate the seroprevalence at each time point.
  sim <- group_by(sim, SAMPLE_ID, ISLAND_ID, TIME) %>%
    summarise(SERO_PREV = sum(R) / sum(S + E + I + R))
  
  # Save the scenario ID and store in the overall data frame
  # of summary metrics.
  all_sims <- sim %>%
    mutate(SCENARIO_ID = i) %>%
    bind_rows(all_sims, .)
}

# Extract the reference values (first scenario).
ref_values <- all_sims %>%
  filter(SCENARIO_ID == 1) %>%
  rename(REF_SERO = SERO_PREV) %>%
  select(-SCENARIO_ID)

# Join in the reference values.
rel_sims <- all_sims %>%
  left_join(ref_values)

# Discard the reference scenario and calculate the difference with reference values.
rel_sims <- ref_sims %>%
  filter(SCENARIO_ID != 1) %>%
  mutate(DIFF_SERO = SERO_PREV - REF_SERO)

# Join with the corresponding dates.
ref_sims <- left_join(ref_sims, dates, by=c("TIME"))

# Restrict dates as requested.
ref_sims <- filter(ref_sims, DATE >= first_date, DATE <= last_date)

# Calculate the median, 50% credible interval and 95% credible interval
# of seroprevalence and infections at each time point for each island.
summary_sims <- ref_sims %>%
  group_by(SCENARIO_ID, ISLAND_ID, DATE) %>%
  summarise(R_500 = quantile(SERO_PREV, probs=0.5),
            R_250 = quantile(SERO_PREV, probs=0.25),
            R_750 = quantile(SERO_PREV, probs=0.75),
            R_975 = quantile(SERO_PREV, probs=0.975),
            R_025 = quantile(SERO_PREV, probs=0.025),
            DIFF_R_500 = quantile(DIFF_SERO, probs=0.5),
            DIFF_R_250 = quantile(DIFF_SERO, probs=0.25),
            DIFF_R_750 = quantile(DIFF_SERO, probs=0.75),
            DIFF_R_975 = quantile(DIFF_SERO, probs=0.975),
            DIFF_R_025 = quantile(DIFF_SERO, probs=0.025))

# Do the same for the reference scenario.
summary_ref <- all_sims %>%
  filter(SCENARIO_ID == 1) %>%
  left_join(dates, by=c("TIME")) %>%
  filter(DATE >= first_date, DATE <= last_date) %>%
  group_by(SCENARIO_ID, ISLAND_ID, DATE) %>%
  summarise(R_500 = quantile(SERO_PREV, probs=0.5),
            R_250 = quantile(SERO_PREV, probs=0.25),
            R_750 = quantile(SERO_PREV, probs=0.75),
            R_975 = quantile(SERO_PREV, probs=0.975),
            R_025 = quantile(SERO_PREV, probs=0.025))

# Define the names of each island based on their ID number.
island_labels <- c("0" = "Anjouan",
                  "1" = "Grande Comore",
                  "2" = "Mayotte",
                  "3" = "Mohéli")

# Define the facet labels.
facet_labels <- names(scenarios)
names(facet_labels) <- paste(seq(1, length(scenarios)))

# Change ISLAND_ID and SCENARIO_ID to a factor.
summary_sims$SCENARIO_ID <- as.factor(summary_sims$SCENARIO_ID)
summary_sims$ISLAND_ID <- factor(summary_sims$ISLAND_ID, levels=c(1, 3, 0, 2))
summary_ref$SCENARIO_ID <- as.factor(summary_ref$SCENARIO_ID)
summary_ref$ISLAND_ID <- factor(summary_ref$ISLAND_ID, levels=c(1, 3, 0, 2))

# Plot seroprevalence over time of the full trade network.
gg_baseline <- ggplot(summary_ref, aes(x=DATE)) +
  geom_ribbon(aes(ymin=R_025, ymax=R_975, fill=ISLAND_ID, alpha="95_CI")) +
  geom_line(aes(y=R_500, colour=ISLAND_ID, alpha="MED"), size=1) +
  scale_x_date(date_breaks="1 year", date_labels="%Y") +
  scale_y_continuous(limits=c(0,1), labels=scales::percent) +
  scale_colour_brewer("Island", palette="Set1",
                      labels=island_labels[levels(summary_ref$ISLAND_ID)]) +
  scale_fill_brewer("Island", palette="Set1",
                    labels=island_labels) +
  scale_alpha_manual("Posterior distribution",
                     breaks=c("MED", "95_CI"),
                     labels=c("MED" = "Median", "95_CI" = "95% credible interval"),
                     values=c("MED" = 0.75, "95_CI" = 0.2)) +
  xlab("Calendar year") +
  ylab("Percentage of animals IgG seropositive") +
  facet_wrap(~SCENARIO_ID, ncol=1, labeller=labeller(SCENARIO_ID=facet_labels)) +
  gg_theme +
  theme(panel.border=element_rect(colour="black", fill=NA)) +
  guides(colour=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                             override.aes=list(alpha=0.5)),
         alpha=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                            override.aes=list(linetype=c(1,0),
                                              alpha=c(1, 0.2),
                                              fill=c(NA, "grey25"))))
print(gg_baseline)

# Plot difference in seroprevalence over time.
gg_diff <- ggplot(summary_sims, aes(x=DATE)) +
  geom_ribbon(aes(ymin=DIFF_R_025, ymax=DIFF_R_975, fill=ISLAND_ID, alpha="95_CI")) +
  geom_line(aes(y=DIFF_R_500, colour=ISLAND_ID, alpha="MED"), size=1) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_date(date_breaks="1 year", date_labels="%Y") +
  scale_y_continuous(limits=c(-1,1), labels=scales::percent) +
  scale_colour_brewer("Island", palette="Set1",
                      labels=island_labels[levels(summary_sims$ISLAND_ID)], guide=FALSE) +
  scale_fill_brewer("Island", palette="Set1",
                    labels=island_labels, guide=FALSE) +
  scale_alpha_manual("Posterior distribution",
                     breaks=c("MED", "95_CI"),
                     labels=c("MED" = "Median", "95_CI" = "95% credible interval"),
                     values=c("MED" = 0.75, "95_CI" = 0.2), guide=FALSE) +
  xlab("Calendar year") +
  ylab("Absolute difference in seroprevalence (%) compared to full trade network") +
  facet_wrap(~SCENARIO_ID, ncol=1, labeller=labeller(SCENARIO_ID=facet_labels)) +
  gg_theme +
  theme(panel.border=element_rect(colour="black", fill=NA))

# Combine the figures.
gg_sp <- egg::ggarrange(gg_baseline, gg_diff, ncol=1, heights=c(1,4.5))

# Save figure as an svg and pdf.
file_name <- "scenario_movement_reduction_seroprev"
ggsave(filename=paste0(file_name, ".svg"), plot=gg_sp, width=16, height=20, dpi=600, scale=0.75,
       path="figures")
ggsave(filename=paste0(file_name, ".pdf"), plot=gg_sp, width=16, height=20, dpi=600, scale=0.75,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, ".pdf figures/", file_name, ".pdf\""))
