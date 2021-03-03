# plot_posterior.R: Visualise the posteriors of the model fit.
# This script takes the posteriors directly from a model fit.

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

# Define the folder where the MCMC output is stored.
mcmc_folder <- "../cpp/out/"

# Set the burn in period.
burn_in <- 600000

# Define the factor to thin the chains by.
thin <- 50

# Define the sample size.
sample_size <- 10000

# Define the number of chains.
n_chains <- 12

# Create space to store all chains.
chains <- data.frame()

# For each chain...
for (i in seq(1, n_chains))
{
  # Load in the MCMC data.
  chain <- fread(paste0(mcmc_folder, "mcmc_", i-1, ".csv"), header=TRUE, sep=",")
  
  # Record the log-likelihood as a parameter.
  log_likelihood <- data.frame("STEP"=chain$STEP[seq(1, nrow(chain), by=length(unique(chain$PAR_NAME)))],
                               "PAR_NAME"="log_likelihood",
                               "PAR_VALUE"=chain$LOG_LIKELIHOOD[seq(1, nrow(chain), by=length(unique(chain$PAR_NAME)))])
  
  # Remove the log-likelihood and prior columns from the chain data.
  chain <- select(chain, -LOG_LIKELIHOOD, -LOG_PRIOR)
  
  # Append this to the chain data frame.
  chain <- rbind(chain, log_likelihood)
  
  # Record the chain ID.
  chain$ID <- i
  
  # Store with all other chains.
  chains <- rbind(chains, chain)
}

# Convert the chain ID to a factor.
chains$ID <- as.factor(chains$ID)

# Get the posterior samples from each chain.
post <- filter(chains, STEP >= burn_in) %>%
  spread(PAR_NAME, PAR_VALUE) %>%
  group_by(ID) %>%
  sample_n(size = sample_size, replace=TRUE) %>%
  gather(key=PAR_NAME, value=PAR_VALUE, -STEP, -ID)

# Thin the chains for plotting.
chains <- chains %>%
  filter(STEP %% thin == 0)

# Load in the prior information.
priors <- fread(paste0(mcmc_folder, "mcmc_priors.csv"), header=TRUE, sep=",")

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
                  "import_start" = "Start of importation window",
                  "import_duration" = "Duration of importation window",
                  "p_imp_inf" = "Proportion of imports that are infectious",
                  "log_likelihood" = "Log-likelihood of data given parameters")

# Construct the trace plots.
gg_chains <- ggplot(chains, aes(x=STEP, y=PAR_VALUE, colour=ID, group=ID)) +
  geom_line(alpha=1/2) +
  scale_colour_brewer("Chain ID", palette="Set1") +
  gg_theme +
  facet_wrap(~PAR_NAME, scales="free_y") +
  theme(panel.border=element_rect(fill=NA, colour="black")) +
  theme(panel.spacing.y = unit(5, "mm")) +
  facet_wrap(~PAR_NAME, scale="free_y", labeller=labeller(PAR_NAME=facet_labels), nrow=5) +
  guides(colour=guide_legend(title.position="top", title.hjust=0, override.aes=list(alpha=1, size=1)))

# Draw the figure.
#plot(gg_chains)

# Construct the posteriors.
gg_post <- post %>%
  ggplot(aes(x = PAR_VALUE, y=..ncount..)) +
  geom_histogram(alpha=0.25, bins=50, col="black", position="identity") +
  stat_density(geom="line", data=priors, aes(x=PAR_VALUE, y=..scaled.., linetype="PRIOR"), adjust=2, col="black") +
  #scale_fill_brewer("Posterior ID", palette="Set1") +
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
file_name <- paste0("model_posterior.svg")
ggsave(filename=file_name, plot=gg_post, width=16, height=9, dpi=600, scale=1,
       path="figures")
system(paste0("bash -c \"pdfcrop --margin '14.22638' figures/", file_name, " figures/", file_name, "\""))

# Get posteriors estimates.
post_est <- post %>%
  filter(PAR_NAME != "log_likelihood") %>%
  group_by(PAR_NAME) %>%
  mutate(PAR_VALUE = ifelse(grepl("move", PAR_NAME), 48*PAR_VALUE, PAR_VALUE),
         PAR_VALUE = ifelse(grepl("import_size", PAR_NAME), 48*PAR_VALUE, PAR_VALUE),
         PAR_VALUE = ifelse(grepl("p_immune", PAR_NAME), 100*PAR_VALUE, PAR_VALUE)) %>%
  summarise(Q_500 = quantile(PAR_VALUE, probs=0.5),
            Q_025 = quantile(PAR_VALUE, probs=0.025),
            Q_975 = quantile(PAR_VALUE, probs=0.975))

# Calculate deviance information criterion (DIC).
# Get the log-liklihood of the best chain.
log_like <- filter(post, PAR_NAME=="log_likelihood")

# Calculate the posterior mean of each parameter.
post_mean <- post %>%
  filter(PAR_NAME != "log_likelihood") %>%
  group_by(PAR_NAME) %>%
  mutate(mu = mean(PAR_VALUE)) %>%
  mutate(diff_mu_par = PAR_VALUE - mu) %>%
  group_by(STEP) %>%
  summarise(DIFF_MU = sum(diff_mu_par^2))

# Find which posterior sample minimises the mean.
mu_post <- post_mean$STEP[which.min(post_mean$DIFF_MU)]

# Get the log-likelihood for the posterior sample which minimises the mean.
log_like_mu_post <- post %>%
  filter(PAR_NAME == "log_likelihood" & STEP==mu_post)
log_like_mu_post <- log_like_mu_post$PAR_VALUE

# Calculate the parameter penalisation factor.
pDICalt <- 2*var(log_like$PAR_VALUE) #alternative.
pDIC <- 2*(log_like_mu_post - mean(log_like$PAR_VALUE))

# Calculate the DIC.
DICalt <- -2*log_like_mu_post + 2*pDICalt
DIC <- -2*log_like_mu_post + 2*pDIC
