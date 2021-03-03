// config.h: Configures macros for directories to read files in from.

#ifndef CONFIG_H
#define CONFIG_H

#define RVF_IDIR "../../data/"  // Relative path (from executable) to the data.
#define RVF_NMCMC 1000000       // Number of MCMC iterations to complete.
#define RVF_CHAINS 12           // How many MCMC chains to run.
#define RVF_INITIAL_PARS 0      // Initialise inferred parameters from pre-selected list (in mcmc.cpp).
#define RVF_ODIR "../out/"      // Relative path (from executable) to the output folder.
#define RVF_RUN_EXTRA 0         // Run the remainder of the time steps after fitting.
#define RVF_SAVE_SIM_FREQ 100   // Frequency to save one simulation trajectory of an accepted parameter set.
#define RVF_UPDATE_FREQ 100     // Frequency to output MCMC information to console and save chains.
#define RVF_UPDATE_COV 1000     // Frequency to update the covariance matrix.
#define RVF_COV_SIZE 10000      // How many elements in the chain to use when calculating mean and covariance.
#define RVF_BEGIN_ADAPT 25000   // Iteration number where adapation begins.
#define RVF_BEGIN_MIX 400000    // Iteration number where chains can mix with one another.
#define RVF_MIX_PROP 0.1        // Proportion of parameter proposals which come from other chains.
#define RVF_MIN_NDVI 0          // Minimum NDVI is local (0) or global (1) or nothing (2).
#define RVF_TRANSMISSION 2      // Use a constant (0), linear (1), or exponential (2) transmission model.
#define RVF_DIFF_SCALE 1        // Use the same (0) or different (1) transmission scalar for islands.
#define RVF_DIFF_NDVI 0         // Use the same (0) or different (1) NDVI scalars for each island.
#define RVF_BIG_NUM 1e10        // Large number to use in place of infinity.

// Define a function which outputs the configuration to file.
void WriteConfig();

#endif