 // likelihood.h: Defines the class Likelihood used to calculate the
 // likelihood of the data given the simulation parameters.
#include "gsl/gsl_rng.h"    // Random number generator.
#include "data.h"           // Definition of data classes.
#include "simulation.h"     // Definition of simulation class.

#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

class Likelihood
{
public:
    // Calculate the log-likelihood.of the data given simulation parameters.
    double LogLikelihood(Data* data, Parameters* pars, Simulation* sim);  

private:
    // Calculate the log-likelihood of a single data entry.
    double SingleLogLikelihood(Data* data, Parameters* pars, Simulation* sim, int entry);
};

#endif