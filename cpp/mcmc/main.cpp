// main.cpp: Describes the workflow of running the RVF (Metras, 2017) model.
#include <iostream>         // Input-output to console.
#include "config.h"         // Number of MCMC chains to run.
#include "data.h"           // Data class definition.
#include "mcmc.h"           // MCMC class definition.
#include "parameters.h"     // Parameter class definition.
#include "likelihood.h"     // Particle filter class definition.
#include "simulation.h"     // Function to simulate RVF dynamics.
#include "omp.h"

int main()
{
    // Output the configuration of the model to file.
    WriteConfig();

    // Space to store parameters for each chain.
    Parameters* pars[RVF_CHAINS];

    // Space to store the data used for the model.
    Data* data[RVF_CHAINS];

    // Space to run the simulations on each chain.
    Simulation* sim[RVF_CHAINS];

    // Space to calculate the likelihood for each chain.
    Likelihood* likelihood[RVF_CHAINS];
    
    // Space to store the MCMC chains.
    Mcmc* mcmc[RVF_CHAINS];

    // Parallelise the entire program to have more than one copy running.
#pragma omp parallel for
    for (int chain_id = 0; chain_id < RVF_CHAINS; ++chain_id)
    {
        // Set up parameters for the model.
        pars[chain_id] = new Parameters{chain_id};
        
        // Set up the data for the model.
        data[chain_id] = new Data{pars[chain_id]};

        // Set up space to run the simulation.
        sim[chain_id] = new Simulation{pars[chain_id], chain_id};

        // Set up space to calculate the likelihood.
        likelihood[chain_id] = new Likelihood{};

        // Set up the MCMC.
        mcmc[chain_id] = new Mcmc{chain_id};

        // Run the MCMC with the current simulation setup.
        // Do a naughty thing and pass the address of all chains, so other
        // chains can be accessed when proposing new parameters.
        mcmc[chain_id]->Run(data[chain_id], pars[chain_id], likelihood[chain_id], sim[chain_id], mcmc);
    }

    // Return an error-free code.
    return 0;
}