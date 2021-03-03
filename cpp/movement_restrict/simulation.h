// simulation.h: Function delcaration for simulating the RVF model.
#include <fstream>          // File stream for recording data.
#include "gsl/gsl_rng.h"    // Random number generation.
#include "data.h"           // Definition of Data class.
#include "parameters.h"     // Definition of Parameters class.

#ifndef SIMULATION_H
#define SIMULATION_H

class Simulation
{
public:
    // Declare variables for the livestock population, split by infectious status:
    // Indexing: time, island and age group.
    double*** S;     // S - Susceptible,
    double*** E;     // E - Exposed,
    double*** I;     // I - Infectious,
    double*** R;     // R - Recovered (immune for life).

    // Some summary measures to keep track of.
    double*** exports;  // exports - Exports of exposed and infectious.
    double*** imports;  // imports - Imports of exposed and infectious into each island over time.

    // Declare a file stream for recording simulation output.
    std::fstream file;

    // Function to simulate the dynamics, including intiialisation,
    // between the desired time points.
    void Simulate(Data* data, Parameters* pars, int start_time, int end_time);

    // Function to write the output of the simulation.
    void WriteOutput(int iter);

    // Constructor and destructor to allocate and de-allocate space for model
    // simulation.
    Simulation(Parameters* pars);
    ~Simulation();

private:
    int n_age;          // Store copy of the number of age groups.
    int n_age_move;     // Store copy of the maximum age group that can be moved.
    int n_island;       // Store copy of the number of islands.
    int n_steps;        // Store copy of the number of time-steps.
    double imp_scaling; // Scaling factor for working out the proportion of imports each age group.
    double*** move_S;   // The number of susceptible movements between islands.
    double*** move_E;   // The number of exposed movements between islands.
    double*** move_I;   // The number of infectious movements between islands.
    double*** move_R;   // The number of recovered movements between islands.

    // Function to calculate the number of movements between islands at
    // a given time step.
    void CalculateMove(Data* data,
                       Parameters* pars,
                       const int t);
};

#endif