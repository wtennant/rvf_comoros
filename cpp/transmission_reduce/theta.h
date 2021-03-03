// theta.h: Defines classes of a parameter to be inferred with different priors.
#include <string>           // String standard library.
#include "gsl/gsl_rng.h"    // Defines random number generator.

#ifndef THETA_H
#define THETA_H

// Base class of an inferred parameter.
class Theta
{
public:
    std::string par_name;               // Name of the parameter being inferred;
    double* chain;                      // Chain of parameter values.
    virtual void Sample(int iter);      // Sample from the underlying parameter distribution.
    Theta(std::string in_par_name);     // Default constructor. Integer to increment seed.
    ~Theta();                           // Default deconstructor.

protected:
    std::string dist_name;  // Name of the parameter distribution.
    int n_pars;             // Number of parameters that define the distribution.
    double* pars;           // Parameter values of the distribution.
};

// Class of an scenario parameter with an uniform distribution.
class Uniform : public Theta
{
public:
    // Function to sample from a uniform distribution.
    void Sample(int scenario_id);

    // Constructor of an inferred parameter with a uniform distribution.
    Uniform(std::string in_par_name, double in_a, double in_b, int inc_seed);

private:
    // Random number generator for sampling from the distribution.
    gsl_rng* rng;
};

#endif  