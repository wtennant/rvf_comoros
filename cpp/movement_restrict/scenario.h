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
    Uniform scale_move_anj_may{"scale_move_anj_may", 0.5, 0.5, 1};              // Scaling of movement between Anjouan and Mayotte.
    Uniform import_duration_two{"import_duration_two", 0.0, 48.0, 2};           // Length of second importation into Grande Comore.
    Uniform import_size_two{"import_size_two", 1.0, 10.0, 3};                   // Number if infectious imports per week in second importation.
    Uniform import_start_two{"import_start_two", 14*48 - 12, 14*48 + 12, 4};    // Start time of the second importation.

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