// likelihood.cpp: Definition of constructor, deconstructor and member functions
// of the Likelihood class. The Likelihood class is used to calculate the
// likelihood of the data given the simulation parameters.
#include <chrono>               // Time for RNG seeding.
#include <iostream>             // Input output stream.
#include "math.h"               // Maximum, logarithm, exponential.
#include "string.h"             // Memory setting.
#include "gsl/gsl_randist.h"    // Resampling particles.
#include "gsl/gsl_sf_gamma.h"   // Evaluating log-likelihood.
#include "config.h"             // Whether to run simulation beyond the data.
#include "likelihood.h"         // Likelihood class definition.

// Function that evaluates the log-likelihood of the data given the simulation
// and its parameters.
double Likelihood::LogLikelihood(Data* data, Parameters* pars, Simulation* sim)
{
    // Initialise the log-likelihood.
    double log_likelihood = 0.0;

    // Initialise the start and end simulation times.
    int start_sim_time = 0;
    int end_sim_time = 0;

    // For each observation...
    for (int i = 0; i < data->sero.n_entries; ++i)
    {
        // Calculate the end simulation time to run the simulation until.
        // If the data entry is aggregated by year...
        if (data->sero.epi_month[i] == 13)
        {
            // ...the final simulation time is the start of the following epi-year.
            end_sim_time = (data->sero.epi_year_start[i] - data->start_year + 1)*12 - data->start_month + 7;
            end_sim_time *= 4;
        }
        else
        {
            // Otherwise, the final simulation time is the start of the following epi-month.
            end_sim_time = (data->sero.epi_year_start[i] - data->start_year)*12 - data->start_month + 7;
            end_sim_time += data->sero.epi_month[i];
            end_sim_time *= 4;            
        }
        
        // Add in a safe guard to ensure the simulation is not run past the
        // available NDVI data.
        end_sim_time = end_sim_time > pars->n_steps ? pars->n_steps : end_sim_time;

        // Run the simulation on all particles between the desired time points.
        if (end_sim_time > start_sim_time)
        {
            sim->Simulate(data, pars, start_sim_time, end_sim_time);
        }

        // Calculate the log-likelihood of the current data entry given the simulation output.
        log_likelihood += SingleLogLikelihood(data, pars, sim, i);
    
        // Start the next simulation at the end of the previous time point.
        start_sim_time = end_sim_time;
    }

    // Should the remainder of the time steps be run?
    if (RVF_RUN_EXTRA == 1)
    {
        // Set the final simulation time to the end time point.
        end_sim_time = pars->n_steps;

        // Only run the remainder of the simulation if the end time exceeds the
        // start time (otherwise the request is pointless!).
        if (end_sim_time > start_sim_time)
        {
            sim->Simulate(data, pars, start_sim_time, end_sim_time);
        }
    }

    // Return the log-likelihood.
    return log_likelihood;
}

// Log-PDF of a binomial distribution with k successes, on n trials with prob
// of success given by p.
double binomial_lpdf(const unsigned int k,
                     const double p, 
                     const unsigned int n)
{
    if (k > n)
    {
        return -1e10;
    }
    else
    {
        double P;
        if (p == 0) 
        {
            P = (k == 0) ? 0 : -1e10;
        }
        else if (p == 1)
        {
            P = (k == n) ? 0 : -1e10;
        }
        else
        {
            double ln_Cnk = gsl_sf_lnchoose(n, k);
            P = ln_Cnk + k * log(p) + (n - k) * log1p(-p);
        }
        return P;
    }
}

// Evaluate the log-likelihood of the current data entry given simulation output.
double Likelihood::SingleLogLikelihood(Data* data, Parameters* pars, Simulation* sim, int entry)
{
    // Calculate the minimum and maximum age for which the seroprevalence
    // needs to be computed given the current data entry.
    // If the age_group entry is -1, age data needs to be aggregated.
    int min_age = (data->sero.age_group[entry] == -1) ? 1: data->sero.age_group[entry];
    int max_age = (data->sero.age_group[entry] == -1) ? 10: data->sero.age_group[entry];

    // If the age_group entry is -2, then all age groups 2--10 should be used.
    min_age = (data->sero.age_group[entry] == -2)? 2: min_age;
    max_age = (data->sero.age_group[entry] == -2)? 10: max_age;

    // If the age_group entry is -3, then this corresponds only to age group 1.
    min_age = (data->sero.age_group[entry] == -3)? 1: min_age;
    max_age = (data->sero.age_group[entry] == -3)? 1: max_age;

    // Calculate the time points in the simulation to aggregate over.
    int start_sim_time;
    int end_sim_time;

    // If the data entry is aggregated by year...
    if (data->sero.epi_month[entry] == 13)
    {
        // ...the final simulation time is the start of the following epi-year.
        start_sim_time = (data->sero.epi_year_start[entry] - data->start_year)*12 - data->start_month + 7;
        end_sim_time = (data->sero.epi_year_start[entry] - data->start_year + 1)*12 - data->start_month + 7;
        start_sim_time *= 4;
        end_sim_time *= 4;
    }
    else
    {
        // Otherwise, the final simulation time is the start of the following epi-month.
        start_sim_time = (data->sero.epi_year_start[entry] - data->start_year)*12 - data->start_month + 7;
        start_sim_time += data->sero.epi_month[entry] - 1;
        end_sim_time = (data->sero.epi_year_start[entry] - data->start_year)*12 - data->start_month + 7;
        end_sim_time += data->sero.epi_month[entry];
        start_sim_time *= 4;
        end_sim_time *= 4;            
    }

    // Add in a safe guard to ensure the simulation is not run past the
    // available NDVI data.
    // Recall one more piece of information is generated than steps in the simulation.
    end_sim_time = (end_sim_time > pars->n_steps)? (pars->n_steps + 1) : end_sim_time;

    // Get the identity of the island the data is from.
    int island_id = data->sero.island_id[entry];

    // Count up the number of recovered and total number of individuals in
    // the requested age bracket and simulation time.
    double recovered = 0.0;
    double total_pop = 0.0;
    for (int t = start_sim_time; t < end_sim_time; ++t)
    {
        for (int age = min_age - 1; age < max_age; ++age)
        {
            recovered += static_cast<double>(sim->R[t][island_id][age]);
            total_pop += static_cast<double>(sim->S[t][island_id][age] +
                sim->E[t][island_id][age] + sim->I[t][island_id][age] +
                    sim->R[t][island_id][age]);
        }
    }

    // Calculate the proportion over the entire time period and age
    // group that are seropositive.
    double prob_sero_pos = recovered / total_pop;

    // Evaluate the probability of observing the data given the
    // seropositivity from the simulation.
    double single_log_likelihood = binomial_lpdf(static_cast<unsigned int>(data->sero.n_positive[entry]),
                    prob_sero_pos, static_cast<unsigned int>(data->sero.n_tested[entry]));

    // Return the log-likelihood of the current data entry.
    return single_log_likelihood;
}
