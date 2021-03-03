 // MCMC.h: Defines the class MCMC used to define priors, likelihood and 
// store parameters to infer using an adaptive Metropolis-Hastings algorithm.
#include "data.h"           // Definition of data classes.
#include "likelihood.h"     // Definition of the likelihood class.
#include "simulation.h"     // Definition of simulation class.
#include "theta.h"          // Definition of MCMC parameter class.

#ifndef MCMC_H
#define MCMC_H

class Mcmc
{
public:
    Mcmc(int chain_id);     // Constructor.
    ~Mcmc();                // Default deconstructor.

    // Run the MCMC algorithm with the given data, model parameters and simulation.
    void Run(Data* data, Parameters* pars, Likelihood* likelihood, Simulation* sim, Mcmc** mcmc);

    // Mean and covariance matrix of the chain.
    double* mean;
    double** covariance;

private:
    // Number of parameters to infer.
    int n_pars;

    // Parameters to infer:
    // Proportion that are immune at time zero for each island.
    Beta p_immune_start_anj{"p_immune_start_anj", 0};
    Beta p_immune_start_gra{"p_immune_start_gra", 1};
    Beta p_immune_start_may{"p_immune_start_may", 2};     
    Beta p_immune_start_moh{"p_immune_start_moh", 3};
    
    // Weekly movements between each island.
    LeftTruncNorm move_anj_gra{"move_anj_gra", 4};
    LeftTruncNorm move_anj_may{"move_anj_may", 5};
    LeftTruncNorm move_anj_moh{"move_anj_moh", 6};
    LeftTruncNorm move_gra_anj{"move_gra_anj", 7};
    LeftTruncNorm move_gra_moh{"move_gra_moh", 8};
    LeftTruncNorm move_moh_anj{"move_moh_anj", 9};
    LeftTruncNorm move_moh_gra{"move_moh_gra", 10};
    
    // Importing parameters.
    LeftTruncNorm import_size{"import_size", 11};
    TruncNorm import_start{"import_start", 12};
    TruncNorm import_duration{"import_duration", 13};

    // Scalar for the NDVI component the transmission rate on each island.
    // Which get used depend on RVF_DIFF_NDVI.
    LeftTruncNorm ndvi_rate{"ndvi_rate", 14};
    LeftTruncNorm ndvi_rate_anj{"ndvi_rate_anj", 15};
    LeftTruncNorm ndvi_rate_gra{"ndvi_rate_gra", 16};
    LeftTruncNorm ndvi_rate_may{"ndvi_rate_may", 17};
    LeftTruncNorm ndvi_rate_moh{"ndvi_rate_moh", 18};

    // Scalar for the constant component of the transmission rate all islands.
    // Which get used depend on RVF_DIFF_SCALE.
    TruncNorm trans_scale{"trans_scale", 19};
    TruncNorm trans_scale_anj{"trans_scale_anj", 20};
    TruncNorm trans_scale_gra{"trans_scale_gra", 21};
    TruncNorm trans_scale_may{"trans_scale_may", 22};
    TruncNorm trans_scale_moh{"trans_scale_moh", 23};

    // Vector used to point to each of the parameters.
    Theta** theta;

    // ID number of the chain.
    int id;

    // Delcare a random number generator.
    gsl_rng* rng;

    // Log-likelihood of the data given the simulation at each step in the chain.
    // Log-prior of the parameters at each step in the chain.
    double* log_likelihood_chain;
    double* log_prior_chain;

    // Acceptance rate of the chain.
    double acceptance_rate;

    // Functions to carry out the MCMC.
    void SetPriors();                               // Sets the priors of the parameters of interest.
    void WritePriors();                             // Output priors to file.
    void Propose(int iter, Mcmc** mcmc);            // Propose a new distribution of parameters.
    void ConstructMeanCov(int iter);                // Calculate the mean and covariance matrix.
    void SetSimPars(int iter, Parameters* pars);    // Set the parameters of the simulation.
    bool IsProposalValid(int iter);                 // Check the validity of proposed parameters.
    double LogPrior(int iter);                      // Evaluate the priors at a given iteration.
    void PrintProgress(int iter);                   // Print progress of the MCMC to console.
    void WriteChains(int iter);                     // Write output of the MCMC.
};

#endif