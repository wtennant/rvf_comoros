# plot_r0.R: Plot the seasonal reproduction number and effective seasonal reproduction number.
# This script uses posteriors and simulations directly from a model fit.

# Clear the workspace.
rm(list = ls())

# Load in data manipulation and visualisation libraries.
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(gridExtra)

# Load in the default ggplot theme.
source("gg_theme.R")

# Define the burn in period of the chains.
burn_in <- 500000

# Load in the MCMC chains.
chain <- fread("../cpp/out/mcmc_1.csv", header=TRUE, sep=",")

# Load in the simulated data.
sim <- fread("../cpp/out/simulation_1.csv", header=TRUE, sep=",")

# Load in the NDVI data.
ndvi <- fread("../data/ndvi.csv", header=TRUE, sep=",") %>%
  mutate(ISLAND_ID = as.numeric(ISLAND_ID)) %>%
  mutate(DATE = as.Date(DATE)) %>%
  group_by(ISLAND_ID) %>%
  mutate(MIN_NDVI = min(NDVI)) %>%
  ungroup()

# Get rid of the burn in period simulations.
sim <- filter(sim, ITER >= burn_in)

# Calculate the proportion susceptible at each time point.
sim <- mutate(sim, WEEK = TIME %% 4,
              MONTH = ((ndvi$MONTH[1] - 1 + (TIME %/% 4)) %% 12) + 1,
              YEAR = ndvi$YEAR[1] + ((ndvi$MONTH[1] - 1 + TIME %/% 4) %/% 12)) %>%
  group_by(ITER, ISLAND_ID, WEEK, MONTH, YEAR) %>%
  summarise(PROP_S = sum(S) / sum(S + E + I + R),
            SERO_PREV = sum(R) / sum(S + E + I + R),
            INFECTIONS = sum(I)) %>%
  rename(STEP=ITER)

# Sample parameter sets that we have the simulated data for.
post <- chain %>%
  filter(STEP >= burn_in) %>%
  spread(key=PAR_NAME, value=PAR_VALUE) %>%
  filter(STEP %in% unique(sim$STEP))

# Calculate the seasonal reproduction number.
rst <- dplyr::select(post, STEP, trans_scale_anj, trans_scale_gra, trans_scale_may, trans_scale_moh, ndvi_rate) %>%
  gather(key=ISLAND_ID, value=trans_scale, contains("trans_scale")) %>%
  mutate(ISLAND_ID = ifelse(ISLAND_ID == "trans_scale_anj", 0, ISLAND_ID)) %>%
  mutate(ISLAND_ID = ifelse(ISLAND_ID == "trans_scale_gra", 1, ISLAND_ID)) %>%
  mutate(ISLAND_ID = ifelse(ISLAND_ID == "trans_scale_may", 2, ISLAND_ID)) %>%
  mutate(ISLAND_ID = ifelse(ISLAND_ID == "trans_scale_moh", 3, ISLAND_ID)) %>%
  mutate(ISLAND_ID = as.numeric(ISLAND_ID)) %>%
  left_join(ndvi, by="ISLAND_ID") %>%
  left_join(sim, by=c("STEP", "ISLAND_ID", "WEEK", "MONTH", "YEAR")) %>%
  mutate(RST = exp(ndvi_rate*(NDVI - MIN_NDVI) + trans_scale))

# Calculate the effective reproduction number.
rest <- rst %>%
  mutate(REST = RST*PROP_S)

# Calculate the median and credible intervals.
sp_ci <- rst %>%
  group_by(ISLAND_ID, DATE) %>%
  summarise(Q_500 = quantile(SERO_PREV, probs=0.5, na.rm=TRUE),
            Q_250 = quantile(SERO_PREV, probs=0.25),
            Q_750 = quantile(SERO_PREV, probs=0.75),
            Q_975 = quantile(SERO_PREV, probs=0.975),
            Q_025 = quantile(SERO_PREV, probs=0.025))
rst_ci <- rst %>%
  group_by(ISLAND_ID, DATE) %>%
  summarise(Q_500 = quantile(RST, probs=0.5),
            Q_250 = quantile(RST, probs=0.25),
            Q_750 = quantile(RST, probs=0.75),
            Q_975 = quantile(RST, probs=0.975),
            Q_025 = quantile(RST, probs=0.025))
rest_ci <- rest %>%
  group_by(ISLAND_ID, DATE) %>%
  summarise(Q_500 = quantile(REST, probs=0.5),
            Q_250 = quantile(REST, probs=0.25),
            Q_750 = quantile(REST, probs=0.75),
            Q_975 = quantile(REST, probs=0.975),
            Q_025 = quantile(REST, probs=0.025))

# Define the names of each island based on their ID number.
facet_labels <- c("0" = "Anjouan",
                  "1" = "Grande Comore",
                  "2" = "Mayotte",
                  "3" = "Mohéli")

# Change the order of islands.
rst$ISLAND_ID <- factor(rst$ISLAND_ID, levels=c(1, 3, 0, 2))
sp_ci$ISLAND_ID <- factor(sp_ci$ISLAND_ID, levels=c(1, 3, 0, 2))
rst_ci$ISLAND_ID <- factor(rst_ci$ISLAND_ID, levels=c(1, 3, 0, 2))
rest_ci$ISLAND_ID <- factor(rest_ci$ISLAND_ID, levels=c(1, 3, 0, 2))

# Plot seroprevalence over time.
gg_sp <- ggplot(sp_ci, aes(x=DATE, colour=ISLAND_ID)) +
  geom_ribbon(aes(ymin=Q_025, ymax=Q_975, fill=ISLAND_ID, alpha="95_CI")) +
  geom_line(aes(y=Q_500, alpha="MED"), size=1) +
  geom_vline(xintercept=as.Date("2015-07-01"), linetype="dashed") +
  #geom_segment(x=as.Date("2015-09-01"), xend=as.Date("2018-01-01"), y=0.95, yend=0.95,
  #             colour="black", arrow=arrow(length=unit(0.5, "cm")), lineend="round", linejoin="mitre") +
  #geom_text(x=as.Date("2015-09-01"), y = 1, vjust=0, hjust=0, label="Forecasting for all islands", colour="black", size=1) +
  geom_vline(xintercept=as.Date("2008-07-01"), linetype="dashed") +
  #geom_segment(x=as.Date("2008-05-01"), xend=as.Date("2006-01-01"), y=0.95, yend=0.95,
  #             colour="black", arrow=arrow(length=unit(0.5, "cm")), lineend="round", linejoin="mitre") +
  #geom_text(x=as.Date("2008-05-01"), y = 1, vjust=0, hjust=1, label="Backcasting for Grande Comore, Mohéli and Anjouan", colour="black", size=1) +
  scale_x_date(date_breaks="1 year", date_labels="%Y") +
  scale_y_continuous(limits=c(0,1), labels=scales::percent) +
  scale_colour_brewer("Island", palette="Set1",
                      labels=facet_labels[levels(rst_ci$ISLAND_ID)]) +
  scale_fill_brewer("Island", palette="Set1",
                    labels=facet_labels) +
  scale_alpha_manual("Posterior distribution",
                     breaks=c("MED", "95_CI"),
                     labels=c("MED" = "Median", "95_CI" = "95% credible interval"),
                     values=c("MED" = 0.75, "95_CI" = 0.2)) +
  xlab("Calendar year") +
  ylab("Percentage of animals seropositive") +
  gg_theme +
  guides(colour=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                             override.aes=list(alpha=0.5)),
         alpha=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                            override.aes=list(linetype=c(1,0),
                                              fill=c(NA, "grey25"))))
print(gg_sp)

# Save and crop.
file_name <- paste0("model_sp_a.pdf")
ggsave(filename=file_name, plot=gg_sp, width=32, height=9, dpi=600, scale=0.6,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, " figures/", file_name, "\""))

# Plot the seasonal reproducction number over time.
gg_rst <- ggplot(rst_ci, aes(x=DATE, colour=ISLAND_ID)) +
  geom_ribbon(aes(ymin=Q_025, ymax=Q_975, fill=ISLAND_ID, alpha="95_CI")) +
  geom_line(aes(y=Q_500, alpha="MED"), size=1) +
  geom_hline(yintercept=1, linetype="dashed", colour="black", alpha=0.75) +
  annotate(geom="text", x=min(rst_ci$DATE), y=1, hjust=1.25, vjust=-0.5, label=expression(paste(R[st]," = 1"))) +
  scale_x_date(date_breaks="1 year", date_labels="%Y") +
  scale_y_continuous(limits=c(0, max(rst_ci$Q_975)), 
                     breaks=seq(0,8,by=1),
                     labels=stringr::str_pad(seq(0,8,by=1), width=4)) +
  scale_colour_brewer("Island", palette="Set1",
                      labels=facet_labels[levels(rst_ci$ISLAND_ID)]) +
  scale_fill_brewer("Island", palette="Set1",
                      labels=facet_labels) +
  scale_alpha_manual("Posterior distribution",
                     breaks=c("MED", "95_CI"),
                     labels=c("MED" = "Median", "95_CI" = "95% credible interval"),
                     values=c("MED" = 0.75, "95_CI" = 0.2)) +
  labs(y=expression(Seasonal~reproduction~number~(R[st]))) +
  xlab("Calendar year") +
  gg_theme +
  theme(legend.text=element_text(colour=NA)) +
  theme(legend.title=element_text(colour=NA)) +
  guides(colour=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                             override.aes=list(alpha=0, colour=NA)),
         alpha=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                            override.aes=list(linetype=c(1,0),
                                              fill=c(NA, "grey50"),
                                              alpha=c(0))))
print(gg_rst)

# Save and crop.
file_name <- paste0("model_rst_a.pdf")
ggsave(filename=file_name, plot=gg_rst, width=32, height=9, dpi=600, scale=0.6,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, " figures/", file_name, "\""))

# Plot the effective seasonal reproducction number over time.
gg_rest <- ggplot(rest_ci, aes(x=DATE, colour=ISLAND_ID)) +
  geom_ribbon(aes(ymin=Q_025, ymax=Q_975, fill=ISLAND_ID, alpha="95_CI")) +
  geom_line(aes(y=Q_500, alpha="MED"), size=1) +
  geom_hline(yintercept=1, linetype="dashed", colour="black", alpha=0.75) +
  annotate(geom="text", x=min(rest_ci$DATE), y=1, hjust=1.25, vjust=-0.5, label=expression(paste(Re[st]," = 1"))) +
  scale_x_date(date_breaks="1 year", date_labels="%Y") +
  scale_y_continuous(limits=c(0, max(rst_ci$Q_975)), breaks=seq(0,8,by=1),
                     labels=stringr::str_pad(seq(0,8,by=1), width=4)) +
  scale_colour_brewer("Island", palette="Set1",
                      labels=facet_labels[levels(rest_ci$ISLAND_ID)]) +
  scale_fill_brewer("Island", palette="Set1",
                    labels=facet_labels) +
  scale_alpha_manual("Posterior distribution",
                     breaks=c("MED", "95_CI"),
                     labels=c("MED" = "Median", "95_CI" = "95% credible interval"),
                     values=c("MED" = 0.75, "95_CI" = 0.2)) +
  labs(y=expression(Effective~seasonal~reproduction~number~(Re[st]))) +
  xlab("Calendar year") +
  gg_theme +
  theme(legend.text=element_text(colour=NA)) +
  theme(legend.title=element_text(colour=NA)) +
  guides(colour=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                             override.aes=list(alpha=0, colour=NA)),
         alpha=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                            override.aes=list(linetype=c(1,0),
                                              fill=c(NA, "grey50"),
                                              alpha=c(0))))
print(gg_rest)

# Combine the plots.
gg_sp_rst_rest <- arrangeGrob(gg_sp, gg_rst, gg_rest)

# Save the figure and crop.
file_name <- paste0("model_rest_a.pdf")
ggsave(filename=file_name, plot=gg_rest, width=32, height=9, dpi=600, scale=0.6,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, " figures/", file_name, "\""))

# Save combined figure and crop.
file_name <- "model_sp_rst_rest"
ggsave(filename=paste0(file_name, ".svg"), plot=gg_sp_rst_rest, width=16, height=9, dpi=600, scale=1.5,
       path="figures")
ggsave(filename=paste0(file_name, ".pdf"), plot=gg_sp_rst_rest, width=16, height=9, dpi=600, scale=1.5,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, ".pdf figures/", file_name, ".pdf\""))

# NDVI against R0.
gg_ndvi_r0 <- group_by(rst, ISLAND_ID, NDVI) %>%
  summarise(Q_500 = quantile(RST, prob=0.5),
            Q_025 = quantile(RST, prob=0.025),
            Q_975 = quantile(RST, prob=0.975)) %>%
  ggplot(aes(x=NDVI, colour=ISLAND_ID)) +
  geom_ribbon(aes(ymin=Q_025, ymax=Q_975, fill=ISLAND_ID, alpha="95_CI")) +
  geom_line(aes(y=Q_500, alpha="MED"), size=1) +
  scale_colour_brewer("Island", palette="Set1",
                      labels=facet_labels[levels(rst$ISLAND_ID)]) +
  scale_fill_brewer("Island", palette="Set1",
                    labels=facet_labels[levels(rst$ISLAND_ID)]) +
  scale_alpha_manual("Posterior predictive distribution",
                     breaks=c("MED", "95_CI"),
                     labels=c("MED" = "Median", "95_CI" = "95% credible interval"),
                     values=c("MED" = 0.75, "95_CI" = 0.2)) +
  xlab("Normalised Difference Vegetation Index (NDVI)") +
  labs(y=expression(Seasonal~reproduction~number~(R[st]))) +
  gg_theme +
  guides(colour=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                             override.aes=list(alpha=0.5)),
         alpha=guide_legend(title.hjust=0, title.position="top", title.vjust=0,
                            override.aes=list(linetype=c(1,0),
                                              fill=c(NA, "grey50"),
                                              alpha=c(1, 0.5))))
print(gg_ndvi_r0)

# Save the figure and crop.
file_name <- paste0("model_ndvi_vs_rst_a.pdf")
ggsave(filename=file_name, plot=gg_ndvi_r0, width=16, height=9, dpi=600, scale=0.75,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, " figures/", file_name, "\""))
