// main.cpp: Describes the workflow of running the RVF (Metras, 2017) model.
#include <iostream>         // Input-output to console.
#include "config.h"         // Number of MCMC chains to run.
#include "data.h"           // Data class definition.
#include "parameters.h"     // Parameter class definition.
#include "posterior.h"      // Posterior class definition.
#include "scenario.h"       // Scenario class definition.
#include "simulation.h"     // Function to simulate RVF dynamics.
#include "omp.h"

int main()
{
    // Output the configuration of the model to file.
    WriteConfig();

    // Set up parameters for the model.
    Parameters pars;
    
    // Set up the data for the model.
    Data data(&pars);

    // Set up space to run the simulation.
    Simulation sim(&pars);
    
    // Set up the posterior distribution to run scenarios from.
    Posterior posterior;

    // Set up the space for scenario testing.
    Scenario scenario;

    // Run the scenario tests.
    scenario.Run(&data, &pars, &posterior, &sim);

    // Return an error-free code.
    return 0;
}