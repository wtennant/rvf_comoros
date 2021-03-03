 // posterior.h: Defines the class Posterior containing all the sampled
 // values from a posterior distribution generated from a model fit.
#include "data.h"           // Definition of data classes.
#include "simulation.h"     // Definition of simulation class.
#include "theta.h"          // Definition of MCMC parameter class.

#ifndef POSTERIOR_H
#define POSTERIOR_H

class Posterior
{
public:
    Posterior();     // Constructor.
    ~Posterior();    // Default deconstructor.

        // Number of parameters in the posterior.
    int n_pars;

    // Number of samples from the posterior.
    int n_samples;

    // All possible parameters in the posterior:
    // Proportion that are immune at time zero for each island.
    Theta p_immune_start_anj{"p_immune_start_anj"};
    Theta p_immune_start_gra{"p_immune_start_gra"};
    Theta p_immune_start_may{"p_immune_start_may"};     
    Theta p_immune_start_moh{"p_immune_start_moh"};
    
    // Weekly movements between each island.
    Theta move_anj_gra{"move_anj_gra"};
    Theta move_anj_may{"move_anj_may"};
    Theta move_anj_moh{"move_anj_moh"};
    Theta move_gra_anj{"move_gra_anj"};
    Theta move_gra_moh{"move_gra_moh"};
    Theta move_moh_anj{"move_moh_anj"};
    Theta move_moh_gra{"move_moh_gra"};
    
    // Importing parameters.
    Theta import_size{"import_size"};
    Theta import_start{"import_start"};
    Theta import_duration{"import_duration"};

    // Scalar for the NDVI component the transmission rate on each island.
    // Which get used depend on RVF_DIFF_NDVI.
    Theta ndvi_rate{"ndvi_rate"};
    Theta ndvi_rate_anj{"ndvi_rate_anj"};
    Theta ndvi_rate_gra{"ndvi_rate_gra"};
    Theta ndvi_rate_may{"ndvi_rate_may"};
    Theta ndvi_rate_moh{"ndvi_rate_moh"};

    // Scalar for the constant component of the transmission rate all islands.
    // Which get used depend on RVF_DIFF_SCALE.
    Theta trans_scale{"trans_scale"};
    Theta trans_scale_anj{"trans_scale_anj"};
    Theta trans_scale_gra{"trans_scale_gra"};
    Theta trans_scale_may{"trans_scale_may"};
    Theta trans_scale_moh{"trans_scale_moh"};

    // Vector used to point to each of the parameters in the posterior.
    Theta** theta;

private:
    // Function to read in and store the posterior distribution.
    void ReadPosterior();

    // Function to write the read in posterior distribution to file (for sanity checking).
    void WritePosterior();
};

#endif