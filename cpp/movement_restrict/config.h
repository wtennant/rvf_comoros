// config.h: Configures macros for directories to read files in from.

#ifndef CONFIG_H
#define CONFIG_H

#define RVF_IDIR "../../data/"  // Relative path (from executable) to the data.
#define RVF_POSTDIR "../in/final/"    // Relative path (from executable) to the folder containing the posterior.
#define RVF_ODIR "out/"         // Relative path (from executable) to the output folder.
#define RVF_UPDATE_FREQ 10      // Frequency to update user on scenario testing progress in console.
#define RVF_NSAMPLES 1000       // Define the number of scenarios to run.
#define RVF_MIN_NDVI 0          // Minimum NDVI is local (0) or global (1) or nothing (2).
#define RVF_TRANSMISSION 2      // Use a constant (0), linear (1), or exponential (2) transmission model.
#define RVF_DIFF_SCALE 1        // Use the same (0) or different (1) transmission scalar for islands.
#define RVF_DIFF_NDVI 0         // Use the same (0) or different (1) NDVI scalars for each island.
#define RVF_SCALE_MAY 0         // Scale movement into Mayotte (1) or not (0).
#define RVF_SEC_IMPORT 0        // Second import into Grande Comore (1) or not (0).
#define RVF_ISOLATE_GRA 1       // Isolate Grande Comore from the network (1) or not (0).
#define RVF_ISOLATE_MOH 0       // Isolate Mohéli from the network (1) or not (0).
#define RVF_ISOLATE_ANJ 0       // Isolate Anjouan from the network (1) or not (0).
#define RVF_ISOLATE_MAY 0       // Isolate Mayotte from the network (1) or not (0).

// Define a function which outputs the configuration to file.
void WriteConfig();

#endif