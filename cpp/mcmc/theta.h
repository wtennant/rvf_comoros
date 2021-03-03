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
    double* chain;                      // Chain containing accepted values of the MCMC.
    virtual double Lpdf(int iter) = 0;  // Virtual function to evaluate the log-probability.
                                        // density function of derived classes.
    virtual bool IsValid(int iter) = 0; // Virtual function to check the validity of a parameter.
    std::string GetInfo();              // Get information of the inferred parameter as a string.
    Theta(int inc_seed);                // Default constructor. Integer to increment seed.
    ~Theta();                           // Default deconstructor.

protected:
    std::string dist_name;  // Name of the prior distribution.
    int n_prior_pars;       // Number of parameters of the prior.
    double* prior_pars;     // Parameter values of the prior.
    gsl_rng* rng;           // Random number generator to sample from the prior.
};

// Class of an inferred parameter with a beta prior.
class Beta : public Theta
{
public:
    // Evaluates log probability of beta distribution.
    double Lpdf(int iter);

    // Check the validity of the parameter at a given iteration in the chain.
    bool IsValid(int iter);

    // Function for setting the prior and generating a first sample.
    void SetPars(double a, double b);

    // Constructor for an inferred parameter with beta prior (a, b).
    Beta(std::string in_par_name, int inc_seed);
};

// Class of an inferred parameter with a normal prior.
class Normal : public Theta
{
public:
    // Evaluates the log probability density function of a normal distribution.
    double Lpdf(int iter);

    // Check the validity of the parameter at a given iteration in the chain.
    bool IsValid(int iter);

    // Function to set the parameters for a prior with a normal distribution.
    void SetPars(double in_mu, double in_sd);

    // Constructor of an inferred parameter with a normal prior.
    Normal(std::string in_par_name, int inc_seed);
};

// Class of an inferrer parameter with a left-truncated normal prior.
class LeftTruncNorm : public Theta
{
public:
    // Evaluates the log probability of a left-truncated normal distribution.
    double Lpdf(int iter);  

    // Check the validity of the parameter at a given iteration in the chain.
    bool IsValid(int iter);

    // Function to set the parameters for a prior with a normal distribution.
    void SetPars(double in_mu, double in_sd, double in_a);

    // Constructor for an inferred parameter with a left-truncated normal distribution
    LeftTruncNorm(std::string in_par_name, int inc_seed);

private:
    double Sample();    // Samples a value from a left-truncated normal distribution.
};

// Class of an inferrer parameter with a right-truncated normal prior.
class RightTruncNorm : public Theta
{
public:
    // Evaluates the log probability of a right-truncated normal distribution.
    double Lpdf(int iter);  

    // Check the validity of the parameter at a given iteration in the chain.
    bool IsValid(int iter);

    // Function to set the parameters for a prior with a right truncated normal distribution.
    void SetPars(double in_mu, double in_sd, double in_b);

    // Constructor for an inferred parameter with a right-truncated normal distribution
    RightTruncNorm(std::string in_par_name, int inc_seed);

private:
    double Sample();    // Samples a value from a left-truncated normal distribution.
};

// Class of an inferrer parameter with a truncated normal prior.
class TruncNorm : public Theta
{
public:
    // Evaluates the log probability of a truncated normal distribution.
    double Lpdf(int iter);  

    // Check the validity of the parameter at a given iteration in the chain.
    bool IsValid(int iter);

    // Function to set the parameters for a prior with a truncated normal distribution.
    void SetPars(double in_mu, double in_sd, double in_a, double in_b);

    // Constructor for an inferred parameter with a left-truncated normal distribution
    TruncNorm(std::string in_par_name, int inc_seed);

private:
    double Sample();    // Samples a value from a left-truncated normal distribution.
};

// Class of an inferred parameter with an uniform prior.
class Uniform : public Theta
{
public:
    // Evaluate the log pdf of a uniform distribution.
    double Lpdf(int iter);

    // Check the validity of the parameter at a given iteration in the chain.
    bool IsValid(int iter);

    // Function to set the parameters for a prior with a uniform distribution.
    void SetPars(double in_a, double in_b);

    // Constructor of an inferred parameter with a uniform prior.
    Uniform(std::string in_par_name, int inc_seed);
};

#endif  