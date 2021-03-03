// theta.cpp: Defines class constructors and member functions of
// inferred parameter classes with a given prior.
#include <chrono>               // Generate random number seed.
#include "math.h"               // Mathematical functions.
#include "config.h"             // Configuration for the length of chains.
#include "gsl/gsl_randist.h"    // Probability density functions.
#include "gsl/gsl_cdf.h"        // Cumulative density functions.
#include "theta.h"              // Class defining parameters to be inferred.
#include <iostream>

// Constructor of the base class Theta.
Theta::Theta(int inc_seed)
{
    // Allocate space to the chain of accepted parameter values.
    chain = new double[RVF_NMCMC + 1]{0.0};
    
    // Set up a random number generator used to sample from the prior.
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    long long seed = std::chrono::duration_cast<std::chrono::nanoseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
    gsl_rng_set(rng, seed + inc_seed);
}

// Default deconstructor of the base class Theta.
Theta::~Theta()
{
    // Free up chain of accepted parameters, parameters for the priors
    // and the random number generator.
    delete[] chain;
    delete[] prior_pars;
    gsl_rng_free(rng);
}

// Function to get information on the inferred parameter as a string.
// This is used to output information on the priors.
std::string Theta::GetInfo()
{
    // Get the parameters as a string, separated by '-'.
    std::string str_prior_pars = "";
    for (int i = 0; i < n_prior_pars; ++i)
    {
        str_prior_pars += (i != 0 ? ";" : "") + std::to_string(prior_pars[i]);
    }

    // Output the concatenation of parameter name, prior distribution name
    // and parameter values used for the prior.
    return par_name + "," + dist_name + "," + str_prior_pars;
}

// Function to evaluates log probability of beta distribution.
double Beta::Lpdf(int iter)
{
    return log(gsl_ran_beta_pdf(chain[iter], prior_pars[0], prior_pars[1]));
}

// Is the parameter at the given iteration in the chain valid?
bool Beta::IsValid(int iter)
{
    return ((chain[iter] >= 0) && (chain[iter] <= 1)) ? true : false;
}

// Function to set the parameters of the beta distribution.
void Beta::SetPars(double in_a, double in_b)
{
    prior_pars[0] = in_a;
    prior_pars[1] = in_b;

    // Randomly sample the first value in the chain of accepted parameters 
    // from the prior.
    chain[0] = gsl_ran_beta(rng, prior_pars[0], prior_pars[1]);
}

// Constructor for the derived class Beta.
Beta::Beta(std::string in_par_name, int inc_seed) : Theta(inc_seed)
{
    // Define the number of parameters in the prior.
    n_prior_pars = 2;

    // Set up space for the parameters of the beta prior distribution.
    prior_pars = new double[n_prior_pars];
    
    // Set up the parameter and prior distribution name.
    par_name = in_par_name; 
    dist_name = "beta";
}

// Function to evaluate the log probability of a normal distribution
// at a given value in the chain of inferred parameter values.
double Normal::Lpdf(int iter)
{
    return log(gsl_ran_gaussian_pdf(chain[iter]  - prior_pars[0], prior_pars[1]));
}

// Is the parameter at the given iteration in the chain valid?
bool Normal::IsValid(int iter)
{
    return true;
}

// Function to set the parameters of the normal distribution.
void Normal::SetPars(double in_mu, double in_sd)
{
    prior_pars[0] = in_mu;
    prior_pars[1] = in_sd;

    // Randomly sample the first value in the chain of accepted parameters 
    // from the prior.
    chain[0] = prior_pars[0] + gsl_ran_gaussian(rng, prior_pars[1]);
}

// Constructor of an inferred parameter with a normal distribution prior.
Normal::Normal(std::string in_par_name, int inc_seed) : Theta(inc_seed)
{
    // Define the number of parameters in the prior.
    n_prior_pars = 2;

    // Set up the parameters of a normal distribution.
    prior_pars = new double[n_prior_pars];

    // Set up the parameter and prior distribution name.
    par_name = in_par_name; 
    dist_name = "norm";
}

// Function to evaluate the log probability density function of a uniform
// distribution at a given value in the chain of inferred parameter values.
double Uniform::Lpdf(int iter)
{
    return log(gsl_ran_flat_pdf(chain[iter], prior_pars[0], prior_pars[1]));
}

// Is the parameter at the given iteration in the chain valid?
bool Uniform::IsValid(int iter)
{
    return ((chain[iter] >= prior_pars[0]) && (chain[iter] <= prior_pars[1])) ? true : false;
}

// Function to set the parameters of the uniform distribution.
void Uniform::SetPars(double in_a, double in_b)
{
    prior_pars[0] = in_a;
    prior_pars[1] = in_b;

    // Randomly sample the first value in the chain of accepted parameters 
    // from the prior.
    chain[0] = gsl_ran_flat(rng, prior_pars[0], prior_pars[1]);
}

// Constructor of an inferred parameter with a uniform prior.
Uniform::Uniform(std::string in_par_name, int inc_seed) : Theta(inc_seed)
{
    // Define the number of parameters in the prior.
    n_prior_pars = 2;

    // Set up the paramters of a uniform distribution.
    prior_pars = new double[n_prior_pars];

    // Set up the parameter and prior distribution name.
    par_name = in_par_name; 
    dist_name = "uniform";
}

// Function to evaluate the log-probability density function
// of the left-truncated normal distribution.
double LeftTruncNorm::Lpdf(int iter)
{
    if (IsValid(iter))
    {
        return log(gsl_ran_gaussian_pdf(chain[iter] - prior_pars[0], prior_pars[1])) -
            log(1.0 - gsl_cdf_gaussian_P(prior_pars[2] - prior_pars[0], prior_pars[1]));
    }
    else
    {
        return -100000.0;   
    }    
}

// Is the parameter at the given iteration in the chain valid?
bool LeftTruncNorm::IsValid(int iter)
{
    return (chain[iter] >= prior_pars[2]) ? true : false;
}

// Function to sample from a truncated normal distribution.
double LeftTruncNorm::Sample()
{
    // Declare the value to be sampled.
    double x;

    // Take the naive strategy of continuously sampling from a Gaussian
    // until inside the bounds.
    do
    {
        x = prior_pars[0] + gsl_ran_gaussian(rng, prior_pars[1]);
    } while(x < prior_pars[2]);

    // Return the sampled value of the truncated normal.
    return x;
}

// Function to set the parameters of the left-truncated distribution.
void LeftTruncNorm::SetPars(double in_mu, double in_sd, double in_a)
{
    prior_pars[0] = in_mu;
    prior_pars[1] = in_sd;
    prior_pars[2] = in_a;

    // Randomly sample the first value in the chain of accepted parameters 
    // from the prior.
    chain[0] = Sample();
}

// Constructor of the left-truncated normal distribution.
LeftTruncNorm::LeftTruncNorm(std::string in_par_name, int inc_seed) : Theta(inc_seed)
{
    // Define the number of parameters in the prior.
    n_prior_pars = 3;

    // Set up the parameters of the priors
    prior_pars = new double[n_prior_pars];

    // Set up the parameter and prior distribution name.
    par_name = in_par_name; 
    dist_name = "ltnorm";
}

// Function to evaluate the log-probability density function
// of the right-truncated normal distribution.
double RightTruncNorm::Lpdf(int iter)
{
    return log(gsl_ran_gaussian_pdf(chain[iter] - prior_pars[0], prior_pars[1])) -
        log(gsl_cdf_gaussian_P(prior_pars[2] - prior_pars[0], prior_pars[1]));
}

// Is the parameter at the given iteration in the chain valid?
bool RightTruncNorm::IsValid(int iter)
{
    return (chain[iter] <= prior_pars[2]) ? true : false;
}

// Function to sample from a truncated normal distribution.
double RightTruncNorm::Sample()
{
    // Declare the value to be sampled.
    double x;

    // Take the naive strategy of continuously sampling from a Gaussian
    // until inside the bounds.
    do
    {
        x = prior_pars[0] + gsl_ran_gaussian(rng, prior_pars[1]);
    } while(x > prior_pars[2]);

    // Return the sampled value of the truncated normal.
    return x;
}

// Function to set the parameters of the right-truncated normal distribution.
void RightTruncNorm::SetPars(double in_mu, double in_sd, double in_b)
{
    prior_pars[0] = in_mu;
    prior_pars[1] = in_sd;
    prior_pars[2] = in_b;

    // Randomly sample the first value in the chain of accepted parameters 
    // from the prior.
    chain[0] = Sample();
}

// Constructor of the right-truncated normal distribution.
RightTruncNorm::RightTruncNorm(std::string in_par_name, int inc_seed) : Theta(inc_seed)
{
    // Define the number of parameters in the prior.
    n_prior_pars = 3;

    // Set up the parameters of the priors
    prior_pars = new double[n_prior_pars];

    // Set up the parameter and prior distribution name.
    par_name = in_par_name; 
    dist_name = "rtnorm";
}

// Function to evaluate the log-probability density function
// of the truncated normal distribution.
double TruncNorm::Lpdf(int iter)
{
    return log(gsl_ran_gaussian_pdf(chain[iter] - prior_pars[0], prior_pars[1])) -
        log(gsl_cdf_gaussian_P(prior_pars[3] - prior_pars[0], prior_pars[1]) -
            gsl_cdf_gaussian_P(prior_pars[2] - prior_pars[0], prior_pars[1]));
}

// Is the parameter at the given iteration in the chain valid?
bool TruncNorm::IsValid(int iter)
{
    return ((chain[iter] <= prior_pars[3]) && (chain[iter] >= prior_pars[2])) ? true : false;
}

// Function to sample from a truncated normal distribution.
double TruncNorm::Sample()
{
    // Declare the value to be sampled.
    double x;

    // Take the naive strategy of continuously sampling from a Gaussian
    // until inside the bounds.
    do
    {
        x = prior_pars[0] + gsl_ran_gaussian(rng, prior_pars[1]);
    } while((x < prior_pars[2]) || (x > prior_pars[3]));

    // Return the sampled value of the truncated normal.
    return x;
}

// Function to set the parameters of the right-truncated normal distribution.
void TruncNorm::SetPars(double in_mu, double in_sd, double in_a, double in_b)
{
    prior_pars[0] = in_mu;
    prior_pars[1] = in_sd;
    prior_pars[2] = in_a;
    prior_pars[3] = in_b;

    // Randomly sample the first value in the chain of accepted parameters 
    // from the prior.
    chain[0] = Sample();
}

// Constructor of the truncated normal distribution.
TruncNorm::TruncNorm(std::string in_par_name, int inc_seed) : Theta(inc_seed)
{
    // Define the number of parameters in the prior.
    n_prior_pars = 4;

    // Set up the parameters of the priors
    prior_pars = new double[n_prior_pars];

    // Set up the parameter and prior distribution name.
    par_name = in_par_name; 
    dist_name = "tnorm";
}