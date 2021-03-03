// theta.cpp: Defines class constructors and member functions of
// the posterior distribution for each parameter.
#include <chrono>               // Generating a seed for a random number generator.
#include <iostream>             // Output to console.
#include "gsl/gsl_randist.h"    // Random number distributions.
#include "theta.h"              // Class defining parameters to be inferred.

// Constructor of the base class Theta.
Theta::Theta(std::string in_par_name)
{
    // Define the parameter name.
    par_name = in_par_name;
}

// Default deconstructor of the base class Theta.
Theta::~Theta()
{
    // Free up chain of posterior parameters.
    delete[] chain;
}

// Sampling from a parameter with no underlying distribution.
void Theta::Sample(int scenario_id)
{
    std::cout << "If you see this message, then something has gone wrong." << std::endl;
    std::cout << "You should not be sampling from a parameter with no underlying distribution!" << std::endl;
}

// Constructor of a parameter with a uniform distribution.
Uniform::Uniform(std::string in_par_name, double in_a, double in_b, int inc_seed) : Theta(in_par_name)
{
    // Define the number of parameters in the distribution.
    n_pars = 2;

    // Set up the paramters of a uniform distribution.
    pars = new double[n_pars]{in_a, in_b};

    // Set up a random number generator used to sample from the distribution.
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    long long seed = std::chrono::duration_cast<std::chrono::nanoseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
    gsl_rng_set(rng, seed + inc_seed);

    // Set up distribution name.
    dist_name = "uniform";
}

// Function to sample from a uniform ditribution and update the list of sampled
// values.
void Uniform::Sample(int scenario_id)
{
    // Randomly sample and save in the chain of parameters.
    chain[scenario_id] = gsl_ran_flat(rng, pars[0], pars[1]);
}