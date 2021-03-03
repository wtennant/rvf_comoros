# plot_posterior_alt.R: Visualise the posteriors of the model fit.
# This script uses a prepared posterior file created by prepare_post.R.

# Clear the workspace.
rm(list = ls())

# Data visulisation and manipulation libraries.
library(dplyr)
library(tidyr)
library(ggplot2)
library(truncnorm)
library(data.table)
library(formattable)
library(gridExtra)
library(tmvtnorm)

# Load in the default ggplot theme.
source("gg_theme.R")

# Define the folder where the prior and posterior is stored.
folder <- "../cpp/in/final/"

# Define size of samples in the prior.
sample_size <- 10000

# Load in the posterior distribution.
post <- fread(paste0(folder, "posterior.csv"), header=TRUE, sep=",")

# Load in the prior information.
priors <- fread(paste0(folder, "mcmc_priors.csv"), header=TRUE, sep=",")

# Separate out the prior parameters.
priors <- suppressWarnings(separate(priors, into=paste0("PAR_", seq(1,4)), col=PRIOR_PARS, sep=";", convert=TRUE))

# Define the priors.
priors <- mutate(priors, PAR_VALUE="0") %>%
  group_by(PAR_NAME) %>%
  mutate(PAR_VALUE = ifelse(PRIOR_DIST=="beta", paste0(rbeta(sample_size, shape1=PAR_1, shape2=PAR_2),collapse=","), PAR_VALUE)) %>%
  mutate(PAR_VALUE = ifelse(PRIOR_DIST=="norm",  paste0(rnorm(sample_size, mean=PAR_1, sd=PAR_2), collapse=","), PAR_VALUE)) %>%
  mutate(PAR_VALUE = ifelse(PRIOR_DIST=="ltnorm", paste0(rtmvnorm(sample_size, mean=PAR_1, sigma=PAR_2^2, lower=PAR_3, upper=Inf), collapse=","), PAR_VALUE)) %>%
  mutate(PAR_VALUE = ifelse(PRIOR_DIST=="rtnorm", paste0(rtmvnorm(sample_size, mean=PAR_1, sigma=PAR_2^2, lower=-Inf, upper=PAR_3), collapse=","), PAR_VALUE)) %>%
  mutate(PAR_VALUE = ifelse(PRIOR_DIST=="tnorm", paste0(rtmvnorm(sample_size, mean=PAR_1, sigma=PAR_2^2, lower=PAR_3, upper=PAR_4), collapse=","), PAR_VALUE)) %>%
  separate(into=paste0("PAR_VALUE_", seq(1,sample_size)), col=PAR_VALUE, sep=",", convert=TRUE) %>%
  gather(key=SAMPLE, value=PAR_VALUE, contains("PAR_VALUE"))

# Define the labels of the facets.
facet_labels <- c("trans_scale" = "Transmission constant",
                  "trans_scale_anj" = "Transmission constant (Anjouan)",
                  "trans_scale_gra" = "Transmission constant (Grande Comore)",
                  "trans_scale_may" = "Transmission constant (Mayotte)",
                  "trans_scale_moh" = "Transmission constant (Mohéli)",
                  "ndvi_rate" = "NDVI transmission coefficient",
                  "ndvi_rate_anj" = "NDVI transmission coefficient (Anjouan)",
                  "ndvi_rate_gra" = "NDVI transmission coefficient (Grande Comore)",
                  "ndvi_rate_may" = "NDVI transmission coefficient (Mayotte)",
                  "ndvi_rate_moh" = "NDVI transmission coefficient (Mohéli)",
                  "p_immune_start_anj" = "Proportion seropositive at time zero (Anjouan)",
                  "p_immune_start_gra" = "Proportion seropositive at time zero (Grande Comore)",
                  "p_immune_start_may" = "Proportion seropositive at time zero (Mayotte)",
                  "p_immune_start_moh" = "Proportion seropositive at time zero (Mohéli)",
                  "move_anj_gra" = "Movement from Anjouan to Grande Comore",
                  "move_anj_may" = "Movement from Anjouan to Mayotte",
                  "move_anj_moh" = "Movement from Anjouan to Mohéli",
                  "move_gra_anj" = "Movement from Grande Comore to Anjouan",
                  "move_gra_moh" = "Movement from Grande Comore to Mohéli",
                  "move_moh_anj" = "Movement from Mohéli to Anjouan",
                  "move_moh_gra" = "Movement from Mohéli to Grande Comore",
                  "import_size" = "Weekly imports into Grande Comore",
                  "import_start" = "Start of importation window (week no. from July 2004)",
                  "import_duration" = "Duration of importation window",
                  "p_imp_inf" = "Proportion of imports that are infectious",
                  "log_likelihood" = "Log-likelihood of data given parameters")

# Construct the posteriors.
gg_post <- post %>%
  ggplot(aes(x = PAR_VALUE, y=..ncount..)) +
  geom_histogram(aes(fill="POSTERIOR"), alpha=0.25, bins=50, col="black", position="identity") +
  stat_density(geom="line", data=priors, aes(x=PAR_VALUE, y=..scaled.., linetype="PRIOR"), adjust=2, col="black") +
  scale_fill_brewer("", palette="Set1", labels="Posterior") +
  scale_linetype_manual("Distribution", labels="Prior", values="dashed") +
  ylab("Scaled density") +
  xlab("Parameter value") +
  gg_theme +
  theme(panel.border=element_rect(fill=NA, colour="black")) +
  theme(legend.title=element_text(margin=margin(0,0,1,0,"mm"))) +
  theme(legend.margin=margin(1,5,1,5,"mm")) +
  theme(panel.spacing.y = unit(5, "mm")) +
  facet_wrap(~PAR_NAME, scale="free_x", labeller=labeller(PAR_NAME=facet_labels), nrow=5) +
  guides(linetype=guide_legend(order=1, title.position="top", title.hjust=0, keywidth=4,
                               override.aes=list(size=1)),
         fill=guide_legend(order=2, title.position="top", title.hjust=0))

# Visualise the posteriors.
plot(gg_post)

# Save the figure and crop.
file_name <- "model_posterior"
ggsave(filename=paste0(file_name, ".svg"), plot=gg_post, width=16, height=12, dpi=600, scale=1.2,
       path="figures")
ggsave(filename=paste0(file_name, ".pdf"), plot=gg_post, width=16, height=12, dpi=600, scale=1.2,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, ".pdf figures/", file_name, ".pdf\""))

# Link the movement parameter names to the island names.
movement_labels <- c("anj" = "Anjouan",
                     "gra" = "Grande Comore",
                     "may" = "Mayotte",
                     "moh" = "Mohéli")

# Get imports and exports on islands from each posterior.
imp_exp <- filter(post, grepl("move", PAR_NAME)) %>%
  separate(col=PAR_NAME, into=c("DISCARD", "SOURCE", "DESTINATION"), sep="_") %>%
  mutate(SOURCE = movement_labels[SOURCE],
         DESTINATION = movement_labels[DESTINATION]) %>%
  select(-DISCARD)

# Calculate total import and exports per island.
exports <- group_by(imp_exp, SOURCE, ID) %>%
  summarise(PAR_VALUE = sum(PAR_VALUE)) %>%
  group_by(SOURCE) %>%
  summarise(Q_500 = 48*quantile(PAR_VALUE, probs=c(0.5)),
            Q_025 = 48*quantile(PAR_VALUE, probs=c(0.025)),
            Q_975 = 48*quantile(PAR_VALUE, probs=c(0.975)))
imports <- group_by(imp_exp, DESTINATION, ID) %>%
  summarise(PAR_VALUE = sum(PAR_VALUE)) %>%
  group_by(DESTINATION) %>%
  summarise(Q_500 = 48*quantile(PAR_VALUE, probs=c(0.5)),
            Q_025 = 48*quantile(PAR_VALUE, probs=c(0.025)),
            Q_975 = 48*quantile(PAR_VALUE, probs=c(0.975)))

# Get some statistics on each edge in the trade network.
move <- filter(post, grepl("move", PAR_NAME)) %>%
  group_by(PAR_NAME) %>%
  summarise(Q_500 = 48*quantile(PAR_VALUE, probs=c(0.5)),
            Q_025 = 48*quantile(PAR_VALUE, probs=c(0.025)),
            Q_975 = 48*quantile(PAR_VALUE, probs=c(0.975)))

# Calculate the posterior median and 95\% credible interval of each parameter.
post_summary <- post %>%
  group_by(PAR_NAME) %>%
  mutate(PAR_VALUE = ifelse(grepl("move", PAR_NAME), 48*PAR_VALUE, PAR_VALUE)) %>%
  mutate(PAR_VALUE = ifelse(grepl("import_size", PAR_NAME), 48*PAR_VALUE, PAR_VALUE)) %>%
  summarise(Q_500 = median(PAR_VALUE),
            Q_025 = quantile(PAR_VALUE, probs=0.025),
            Q_975 = quantile(PAR_VALUE, probs=0.975))
