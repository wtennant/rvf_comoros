// Scenario.h: Defines the class Scenario used to run different scenario tests
// given the posterior from a previous model fit.
#include "data.h"           // Definition of Data class.
#include "posterior.h"      // Definition of the Posterior class.
#include "simulation.h"     // Definition of Simulation class.
#include "theta.h"          // Definition of MCMC parameter class.

#ifndef SCENARIO_H
#define SCENARIO_H

class Scenario
{
public:
    Scenario();                 // Default constructor.
    ~Scenario();                // Default deconstructor.

    // Function to start running the scenario tests.
    void Run(Data* data, Parameters* pars, Posterior* post, Simulation* sim);

private:
    // Number of scenarios to implement simultaneously.
    int n_scenarios;

    // Declare a random number generator.
    gsl_rng* rng;

    // Sample IDs for the posteriors used.
    int* post_idx;

    // Distributions for parameters in scenarios that are to be tested.
    Uniform vec_control_gra{"vec_control_gra", 0.0, 0.0, 1};            // Vector control on Grande Comore.
    Uniform vec_control_moh{"vec_control_moh", 0.1, 0.3, 2};            // Vector control on Mohéli.
    Uniform vec_control_anj{"vec_control_anj", 0.0, 0.0, 3};            // Vector control on Anjouan.
    Uniform vec_control_may{"vec_control_may", 0.0, 0.0, 4};            // Vector control on Mayotte

    // Space to store all scenario parameters.
    Theta** theta;

    // Functions to carry out the scenario tests.
    void SetSimPars(int sample_id, Parameters* pars, Posterior* post);  // Set the parameters of the simulation from the posterior distribution.
    void UpdateSimPars(int sample_id, Parameters* pars);                // Update the simulation parameter sets based on sampled scenarios.
    void PrintProgress(int sample_id);                                  // Print progress of scenario tests to console.
    void WritePosteriorSamples(Posterior* post);                        // Output sampled posterior parameters to file.
    void WriteScenarioPars();                                           // Output scenario parameters to file.
};

#endif