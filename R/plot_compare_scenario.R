# plot_compare_scenario.R: Compares all isolation scenarios.

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
folder <- "../cpp/movement_restrict/out/"

# Define a look up table with the scenario name and the
# sub-folder its contained in
scenarios <- c("Full trade network" = "scenario_free/",
               "Grande Comore isolation" = "isolate_grande_comore/",
               "Mohéli isolation" = "isolate_moheli/",
               "Anjouan isolation" = "isolate_anjouan/",
               "Mayotte isolation" = "isolate_mayotte/")

# Create space to store the summary metrics of all scenarios.
sum_scenario <- data.frame()

# For each scenario, calculate the total number of infections
# on each island for each simulation.
for (i in seq(1, length(scenarios)))
{
  # Load in the simulation data.
  sim <- fread(paste0(folder, scenarios[i], "simulation.csv"), header=TRUE)
  
  # Group by island and simulation ID, and sum up all infections.
  sim <- group_by(sim, SAMPLE_ID, ISLAND_ID) %>%
    filter(TIME < 528) %>%
    summarise(INF = sum(E+I)/2)
  
  # Save the scenario ID and store in the overall data frame
  # of summary metrics.
  sum_scenario <- sim %>%
    mutate(SCENARIO_ID = i) %>%
    bind_rows(sum_scenario, .)
}

# Get the median of all simulations in the scenario-free model.
scen_free <- filter(sum_scenario, SCENARIO_ID == 1) %>%
  rename(SCEN_FREE_INF = INF) %>%
  select(ISLAND_ID, SAMPLE_ID, SCEN_FREE_INF)

# Join onto the scenario summary data frame. 
sum_scenario <- left_join(sum_scenario, scen_free, by=c("ISLAND_ID","SAMPLE_ID"))

# Compute INF relative to the scenario free sim.
sum_scenario <- mutate(sum_scenario, REL_INF = INF - SCEN_FREE_INF)

# Convert scenario and island IDs to factors.
sum_scenario$SCENARIO_ID <- as.factor(sum_scenario$SCENARIO_ID)
sum_scenario$ISLAND_ID <- factor(sum_scenario$ISLAND_ID, levels=c(1, 3, 0, 2))
sum_scenario$FACET_ORDER <- factor(sum_scenario$ISLAND_ID, levels=c(1, 0, 3, 2))

# Calculate infections across all islands.
sum_total_inf <- sum_scenario %>%
  group_by(SAMPLE_ID, SCENARIO_ID) %>%
  summarise(INF = sum(INF))

# Linking ISLAND_IDs to island names.
island_name <- c("1" = "Grande Comore",
                 "3" = "Mohéli",
                 "0" = "Anjouan",
                 "2" = "Mayotte")

# For each scenario, make a boxplot.
gg_scenario <- ggplot(sum_scenario, aes(x=SCENARIO_ID, y=INF)) +
  #geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) +
  geom_violin(aes(fill=ISLAND_ID, group=as.factor(paste(SCENARIO_ID, as.numeric(ISLAND_ID)))), scale="width", position=position_dodge(0.8), alpha=0.75, width = 0.7) +
  geom_boxplot(aes(group=as.factor(paste(SCENARIO_ID, as.numeric(ISLAND_ID)))), alpha=0.75, position=position_dodge(0.8), width=0.3) +
  scale_x_discrete("Movement restriction scenario", labels=names(scenarios)) +
  scale_y_continuous("Number of infections from July 2004 to June 2015", labels=scales::comma) +
  scale_fill_brewer("Island", palette="Set1", labels=c("Grande Comore", "Mohéli", "Anjouan", "Mayotte"), guide=FALSE) +
  gg_theme +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  theme(panel.border=element_rect(colour="black", fill=NA)) +
  facet_wrap(~FACET_ORDER, scales="free_y", labeller=labeller(FACET_ORDER=island_name)) +
  theme(aspect.ratio=1)

# Create a 1x1 panel plot of the total number of infections for each scenario.
gg_total_inf <- ggplot(sum_total_inf, aes(x=SCENARIO_ID, y=INF)) +
  geom_violin(aes(group=SCENARIO_ID), scale="width", position=position_dodge(0.8), alpha=0.75, width = 0.7, fill="grey85") +
  geom_boxplot(aes(group=SCENARIO_ID), alpha=0.75, position=position_dodge(0.8), width=0.3, fill="white") +
  scale_x_discrete("Movement restriction scenario", labels=names(scenarios)) +
  scale_y_continuous("Number of infections from July 2004 to June 2015", labels=scales::comma) +
  gg_theme +
  facet_wrap(~"Comoros archipelago") +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  theme(aspect.ratio=1) +
  theme(panel.border=element_rect(colour="black", fill=NA))

# Combine the figures.
p1 <- egg::ggarrange(gg_total_inf, gg_scenario, ncol=2)
print(p1)

# Save figure as an svg and pdf.
file_name <- "scenario_movement_reduction"
ggsave(filename=paste0(file_name, ".svg"), plot=p1, width=16, height=9, dpi=600, scale=1.0,
       path="figures")
ggsave(filename=paste0(file_name, ".pdf"), plot=p1, width=16, height=9, dpi=600, scale=1.0,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, ".pdf figures/", file_name, ".pdf\""))

# Summarise the summary...
sum_sum_scenario <- sum_scenario %>%
  group_by(SCENARIO_ID, ISLAND_ID) %>%
  summarise(Q_500 = quantile(100*REL_INF/SCEN_FREE_INF, prob=0.5),
            Q_975 = quantile(100*REL_INF/SCEN_FREE_INF, prob=0.975),
            Q_025 = quantile(100*REL_INF/SCEN_FREE_INF, prob=0.025))
sum_inf_scenario <- sum_scenario %>%
  group_by(SCENARIO_ID, SAMPLE_ID) %>%
  summarise(INF = sum(INF)) %>%
  summarise(Q_500 = quantile(INF, prob=0.5),
            Q_975 = quantile(INF, prob=0.975),
            Q_025 = quantile(INF, prob=0.025))
rel_inf_scenario <- sum_scenario %>%
  group_by(SCENARIO_ID, SAMPLE_ID) %>%
  summarise(REL_INF = sum(INF - SCEN_FREE_INF)) %>%
  summarise(Q_500 = quantile(REL_INF, prob=0.5),
            Q_975 = quantile(REL_INF, prob=0.975),
            Q_025 = quantile(REL_INF, prob=0.025))
