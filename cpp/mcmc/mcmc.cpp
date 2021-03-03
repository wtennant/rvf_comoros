// mcmc.cpp: Defines the function to run the Metropolis-Hastings MCMC algorithm.
#include <fstream>              // Writing to file.
#include <iostream>             // Output information to console.
#include <chrono>               // Time used as the random number generator seed.
#include <iomanip>              // Set precision on file output.
#include "math.h"               // Mathematical functions including log and exp.
#include "gsl/gsl_linalg.h"     // Cholesky decomposition for the proposal distribution.
#include "gsl/gsl_randist.h"    // Distributions to evaluate likelihood.
#include "config.h"             // Configuration of file output.
#include "mcmc.h"               // MCMC class definition.

// Function to set the simulation parameters.
// Note: when new parameters (to infer) are added or removed, this must be
// changed. Default parameters are set in parameters.cpp.
void Mcmc::SetSimPars(int iter, Parameters* pars)
{
    // Set the parameters in the simulation to be the
    // proposed set of parameters in theta.
    pars->p_immune_start[0] = p_immune_start_anj.chain[iter];
    pars->p_immune_start[1] = p_immune_start_gra.chain[iter];
    pars->p_immune_start[2] = p_immune_start_may.chain[iter];
    pars->p_immune_start[3] = p_immune_start_moh.chain[iter];
    pars->move[0][1] = move_anj_gra.chain[iter];
    pars->move[0][2] = move_anj_may.chain[iter];
    pars->move[0][3] = move_anj_moh.chain[iter];
    pars->move[1][0] = move_gra_anj.chain[iter];
    pars->move[1][3] = move_gra_moh.chain[iter];
    pars->move[3][0] = move_moh_anj.chain[iter];
    pars->move[3][1] = move_moh_gra.chain[iter];
    pars->import_start = import_start.chain[iter];
    pars->import_duration = import_duration.chain[iter];
    pars->import_size = import_size.chain[iter];

    // The model parameters used depend on configuration.
    switch(RVF_DIFF_NDVI)
    {
        case 0:     // Same NDVI per island.
            pars->ndvi_rate[0] = ndvi_rate.chain[iter];
            pars->ndvi_rate[1] = ndvi_rate.chain[iter];
            pars->ndvi_rate[2] = ndvi_rate.chain[iter];
            pars->ndvi_rate[3] = ndvi_rate.chain[iter];
            break;
        case 1:     // Different NDVI per island.
            pars->ndvi_rate[0] = ndvi_rate_anj.chain[iter];
            pars->ndvi_rate[1] = ndvi_rate_gra.chain[iter];
            pars->ndvi_rate[2] = ndvi_rate_may.chain[iter];
            pars->ndvi_rate[3] = ndvi_rate_moh.chain[iter];
            break;
        default:    // Default is same.
            pars->ndvi_rate[0] = ndvi_rate.chain[iter];
            pars->ndvi_rate[1] = ndvi_rate.chain[iter];
            pars->ndvi_rate[2] = ndvi_rate.chain[iter];
            pars->ndvi_rate[3] = ndvi_rate.chain[iter];
    }
    switch(RVF_DIFF_SCALE)
    {
        case 0:     // Same transmission scale per island.
            pars->trans_scale[0] = trans_scale.chain[iter];
            pars->trans_scale[1] = trans_scale.chain[iter];
            pars->trans_scale[2] = trans_scale.chain[iter];
            pars->trans_scale[3] = trans_scale.chain[iter];
            break;
        case 1:     // Different transmission scale per island.
            pars->trans_scale[0] = trans_scale_anj.chain[iter];
            pars->trans_scale[1] = trans_scale_gra.chain[iter];
            pars->trans_scale[2] = trans_scale_may.chain[iter];
            pars->trans_scale[3] = trans_scale_moh.chain[iter];
            break;
        default:    // Default is same.
            pars->trans_scale[0] = trans_scale.chain[iter];
            pars->trans_scale[1] = trans_scale.chain[iter];
            pars->trans_scale[2] = trans_scale.chain[iter];
            pars->trans_scale[3] = trans_scale.chain[iter];
    }
        
}

// Set all the priors for each parameter.
void Mcmc::SetPriors()
{
    // Initialisation and movement priors.
    p_immune_start_anj.SetPars(5.0, 45.0);                          // Immunity in Anjouan at time zero.
    p_immune_start_gra.SetPars(20.0, 30.0);                         // Immunity in Grande Comore at time zero.
    p_immune_start_may.SetPars(5.0, 45.0);                          // Immunity in Mayotte at time zero.
    p_immune_start_moh.SetPars(20.0, 30.0);                         // Immunity in Mohéli at time zero.
    move_anj_gra.SetPars(625.0 / 48.0, 625.0 / 48.0 / 10.0, 1.0);   // Weekly movement from A to GC.
    move_anj_may.SetPars(1500.0 / 48.0, 1500.0 / 48.0 / 10.0, 1.0); // Weekly movement from A to Ma.
    move_anj_moh.SetPars(409.0 / 48.0, 409.0 / 48.0 / 10.0, 1.0);   // Weekly movement from A to Mo.
    move_gra_anj.SetPars(50.0 / 48.0, 50.0 / 48.0, 0.0);            // Weekly movement from GC to A.
    move_gra_moh.SetPars(350.0 / 48.0, 350.0 / 48.0 / 10.0, 1.0);   // Weekly movement from GC to Mo.
    move_moh_anj.SetPars(50.0 / 48.0, 50.0 / 48.0, 0.0);            // Weekly movement from Mo to A.
    move_moh_gra.SetPars(676.0 / 48.0, 676.0 / 48.0 / 10.0, 1.0);   // Weekly movement from Mo to GC.
    import_size.SetPars(150.0 / 48.0, 150.0 / 48.0, 0.0);           // Weekly infections in import window.
    import_start.SetPars(128.0, 24.0, 116.0, 156.0);                // Start of import window.
    import_duration.SetPars(24.0, 12.0, 0.0, 48.0);                 // Duration of import window.

    // Transmission parameters.
    ndvi_rate.SetPars(3.0, 2.0, 0.0);  
    ndvi_rate_anj.SetPars(3.0, 2.0, 0.0);                           
    ndvi_rate_gra.SetPars(3.0, 2.0, 0.0);
    ndvi_rate_may.SetPars(3.0, 2.0, 0.0);
    ndvi_rate_moh.SetPars(3.0, 2.0, 0.0);
    
    // The boundaries of the constant transmission rate depend on the model.
    switch (RVF_TRANSMISSION)
    {
        case 0:     // Constant transmission model.
            trans_scale.SetPars(1.0, 2.0, 0.0, RVF_BIG_NUM);
            trans_scale_anj.SetPars(1.0, 2.0, 0.0, RVF_BIG_NUM);
            trans_scale_gra.SetPars(1.0, 2.0, 0.0, RVF_BIG_NUM);
            trans_scale_may.SetPars(1.0, 2.0, 0.0, RVF_BIG_NUM);
            trans_scale_moh.SetPars(1.0, 2.0, 0.0, RVF_BIG_NUM);
            break;
        case 1:     // Linear transmission model.
            trans_scale.SetPars(1.0, 2.0, 0.0, RVF_BIG_NUM);
            trans_scale_anj.SetPars(1.0, 2.0, 0.0, RVF_BIG_NUM);
            trans_scale_gra.SetPars(1.0, 2.0, 0.0, RVF_BIG_NUM);
            trans_scale_may.SetPars(1.0, 2.0, 0.0, RVF_BIG_NUM);
            trans_scale_moh.SetPars(1.0, 2.0, 0.0, RVF_BIG_NUM);
            break;
        case 2:     // Exponential transmission model.
        
            trans_scale.SetPars(-2.0, 2.0, -RVF_BIG_NUM, 0.0);
            trans_scale_anj.SetPars(-2.0, 2.0, -RVF_BIG_NUM, 0.0);
            trans_scale_gra.SetPars(-2.0, 2.0, -RVF_BIG_NUM, 0.0);
            trans_scale_may.SetPars(-2.0, 2.0, -RVF_BIG_NUM, 0.0);
            trans_scale_moh.SetPars(-2.0, 2.0, -RVF_BIG_NUM, 0.0);
            break;
        default:    // Default to the exponential.
            trans_scale.SetPars(-2.0, 2.0,  -RVF_BIG_NUM, 0.0);
            trans_scale_anj.SetPars(-2.0, 2.0,  -RVF_BIG_NUM, 0.0);
            trans_scale_gra.SetPars(-2.0, 2.0,  -RVF_BIG_NUM, 0.0);
            trans_scale_may.SetPars(-2.0, 2.0,  -RVF_BIG_NUM, 0.0);
            trans_scale_moh.SetPars(-2.0, 2.0,  -RVF_BIG_NUM, 0.0);
    }
    
}

// Constructor of the Mcmc class.
Mcmc::Mcmc(int chain_id)
{
    // Define the number of parameters.
    n_pars = 14;    // The number of parameters independent of model choice.

    // Add in the number of parameters depending on:
    // different ndvi influence on transmission per island (or not),
    // and different transmission scales per island (or not).
    switch(RVF_DIFF_NDVI)
    {
        case 0: // NDVI influence is the same on each island.
            n_pars += 1;
            break;
        case 1: // NDVI influence is different on each island.
            n_pars += 4;
            break;
        default:
            n_pars += 1;
    }
    switch(RVF_DIFF_SCALE)
    {
        case 0: // Transmission scale is the same on each island.
            n_pars += 1;
            break;
        case 1: // Transmission scale is different on each island.
            n_pars += 4;
            break;
        default:
            n_pars += 1;
    }

    // Initialise the acceptance rate.
    acceptance_rate = 0.0;

    // Store a local copy of the chain id.
    id = chain_id;

    // Create space to store the log-likelihood, log-prior and log-posterior
    // at each iteration in the chain.
    log_likelihood_chain = new double[RVF_NMCMC + 1];
    log_prior_chain = new double[RVF_NMCMC + 1];

    // Allocate space for the empirical mean and covariance.
    mean = new double[n_pars]{0.0};
    covariance = new double*[n_pars];
    for (int i = 0; i < n_pars; ++i)
    {
        covariance[i] = new double[n_pars]{0.0};
    }

    // Set up the random number generator based on time and chain ID.
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, 2*std::chrono::system_clock::now().time_since_epoch().count() + chain_id);

    // Set all the priors.
    SetPriors();

    // Space to store all the parameters and their chains.
    theta = new Theta*[n_pars]{&p_immune_start_anj,
                               &p_immune_start_gra,
                               &p_immune_start_may,
                               &p_immune_start_moh,
                               &move_anj_gra,
                               &move_anj_may,
                               &move_anj_moh,
                               &move_gra_anj,
                               &move_gra_moh,
                               &move_moh_anj,
                               &move_moh_gra,
                               &import_size,
                               &import_start,
                               &import_duration};

    // This will depend on the model: whether NDVI influence, and/or transmission
    // scale, are the same or different on each island.
    switch(RVF_DIFF_NDVI)
    {
        case 1:
            theta[14] = &ndvi_rate_anj;
            theta[15] = &ndvi_rate_gra;
            theta[16] = &ndvi_rate_may;
            theta[17] = &ndvi_rate_moh;
            switch(RVF_DIFF_SCALE)
            {
                case 1:
                    theta[18] = &trans_scale_anj;
                    theta[19] = &trans_scale_gra;
                    theta[20] = &trans_scale_may;
                    theta[21] = &trans_scale_moh;
                    break;
                default:
                    theta[18] = &trans_scale;
            }
            break;
        default:
            theta[14] = &ndvi_rate;
            switch(RVF_DIFF_SCALE)
            {
                case 1:
                    theta[15] = &trans_scale_anj;
                    theta[16] = &trans_scale_gra;
                    theta[17] = &trans_scale_may;
                    theta[18] = &trans_scale_moh;
                    break;
                default:
                    theta[15] = &trans_scale;
            }
    }
    
    // Write the priors to file. Only do this once.
    if (chain_id == 0)
    {
        WritePriors();
    }
}

// Deconstructor of the Mcmc class.
Mcmc::~Mcmc()
{
    // De-allocate the random number generator.
    gsl_rng_free(rng);

    // Free up the parameter and log-likelihood space.
    for (int i = 0; i < n_pars; ++i)
    {
        delete[] covariance[i];
    }
    delete[] covariance;
    delete[] mean;
    delete[] log_likelihood_chain;
    delete[] theta;
}

// Metropolis-Hastings MCMC algorithm.
void Mcmc::Run(Data* data, Parameters* pars, Likelihood* likelihood, Simulation* sim, Mcmc** mcmc)
{
    // Initialise the inferred parameters if requested.
    if (RVF_INITIAL_PARS == 1)
    {
        p_immune_start_anj.chain[0] = 0.0286;
        p_immune_start_gra.chain[0] = 0.386;
        p_immune_start_may.chain[0] = 0.150;
        p_immune_start_moh.chain[0] = 0.463;
        trans_scale_anj.chain[0] = -1.12;
        trans_scale_gra.chain[0] = -0.797;
        trans_scale_may.chain[0] = -1.08;
        trans_scale_moh.chain[0] = -0.838;
        ndvi_rate.chain[0] = 8.16;
        move_anj_gra.chain[0] = 16.1;
        move_anj_may.chain[0] = 25.6;
        move_anj_moh.chain[0] = 6.51;
        move_gra_anj.chain[0] = 1.57;
        move_gra_moh.chain[0] = 9.21;
        move_moh_anj.chain[0] = 4.26;
        move_moh_gra.chain[0] = 14.4;
        import_size.chain[0] = 9.46;
        import_start.chain[0] = 117;
        import_duration.chain[0] = 13.2;
    }

    // Set the simulation parameters on the parameter values sampled
    // from the prior.
    SetSimPars(0, pars);

    // Calculate the log-likelihood of the initial parameter set.
    // Calculate the log-prior of the intiial parameter set.
    log_likelihood_chain[0] = likelihood->LogLikelihood(data, pars, sim);
    log_prior_chain[0] = LogPrior(0);

    // Write output to file.
    sim->WriteOutput(0);

    // Delcare a variable denoting whether or not a proposed parameter set
    // should be rejected.
    bool reject;

    // For each iteration in the MCMC.
    for (int iter = 0; iter < RVF_NMCMC; ++iter)
    {
        // Propose a new parameter set based on the previously accepted values.
        Propose(iter, mcmc);

        // Set the simulation parameters to the parameters sampled from the
        // proposal distribution.
        SetSimPars(iter + 1, pars);

        // Run a particle filter if the proposed set of parameters are valid.
        if (IsProposalValid(iter + 1))
        {
            // Calculate the log-likelihood of the proposed parameters.
            log_likelihood_chain[iter + 1] = likelihood->LogLikelihood(data, pars, sim);

            // Calculate the log-prior of the proposed parameters.
            log_prior_chain[iter + 1] = LogPrior(iter + 1);

            // Calculate the log acceptance ratio of the proposed parameter set.
            // Note that as the proposal distribution is symmetric, there is no need
            // to include the term account for bias in the proposal distribution.
            double log_acceptance_ratio = log_likelihood_chain[iter + 1] - log_likelihood_chain[iter];
            log_acceptance_ratio += log_prior_chain[iter + 1] - log_prior_chain[iter]; // Factor in the priors.

            // Generate a random number between zero and one, and compare with the
            // acceptance ratio. Reject if acceptance ratio is too small.
            reject = (gsl_ran_flat(rng, 0, 1) > exp(log_acceptance_ratio)) ? true : false;
        }
        else
        {
            // Reject the parameter set if the proposal is invalid.
            reject = true;
        }

        // Update the chain depending on whether the proposed parameter set was
        // rejected or accepted.
        if (reject)
        {
            // Reject the proposed parameter set.
            // For each parameter, set the current iteration to the previous.
            for (int i = 0; i < n_pars; ++i)
            {
                theta[i]->chain[iter+ 1] = theta[i]->chain[iter];
            }

            // Do the same for the log-likelihood.
            log_likelihood_chain[iter + 1] = log_likelihood_chain[iter];
            log_prior_chain[iter + 1] = log_prior_chain[iter];

            // Update the acceptance ratio.
            acceptance_rate = acceptance_rate * static_cast<double>(iter) /
                static_cast<double>(iter + 1);
        }
        else
        {
            // Accept the proposed parameter set.
            // By construction, the only thing to do is to
            // update the acceptance ratio.
            acceptance_rate = acceptance_rate * static_cast<double>(iter) /
                static_cast<double>(iter + 1) + 1.0 / static_cast<double>(iter + 1);
        }
        
        // Print an update to console. Only print progress on the master chain.
        if (id == 0)
        {
            PrintProgress(iter + 1);
        }

        // Write chains to file.
        if ((iter + 1) % RVF_UPDATE_FREQ == 0)
        {
            WriteChains(iter + 1);
        }
    }

    // Write the remainder of the chains to file.
    WriteChains(RVF_NMCMC + 1);
}

// Function to calculate the mean and covariance of the chain.
void Mcmc::ConstructMeanCov(int iter)
{
    // Define the first iteration to use in calculating the mean
    // and covariance of each chain.
    double first_iter = iter < RVF_COV_SIZE ? 0 : iter - RVF_COV_SIZE + 1;

    // Construct the mean of each chain.
    for (int i = 0; i < n_pars; ++i)
    {
        mean[i] = 0.0;
        for (int k = first_iter; k <= iter; ++k)
        {
            mean[i] += theta[i]->chain[k] / static_cast<double>(iter - first_iter + 1);
        }

        for (int j = i; j < n_pars; ++j)
        {
            covariance[i][j] = 0.0;
            for (int k = first_iter; k <= iter; ++k)
            {
                covariance[i][j] += (theta[i]->chain[k] - mean[i]) *
                    (theta[j]->chain[k] - mean[j]) / static_cast<double>(iter - first_iter);
            }
            covariance[j][i] = covariance[i][j];
        }
    }
}   

// Function to propose a new parameter set based on an adaptive random-walk.
void Mcmc::Propose(int iter, Mcmc** mcmc)
{
    // Here, we propose parameters based on two distributions.
    // The first will be a random walk from the chain itself. 
    // The second will be proposing from other chains that are running.
    if ((iter < RVF_BEGIN_MIX) || (gsl_ran_flat(rng, 0, 1) > RVF_MIX_PROP) || (RVF_CHAINS == 1))
    {
        // Construct the small random-walk matrix:
        // (0.1)^2 I / n_pars
        gsl_matrix* cov = gsl_matrix_alloc(n_pars, n_pars);
        gsl_matrix_set_identity(cov);
        gsl_matrix_scale(cov, 0.1 * 0.1 / static_cast<double>(n_pars));

        // Set the weight to the empirical covariance matrix versus fixed random walk.
        double beta_prob = 0.05;

        // Use the classic Roberts and Rosenthal definition for a proposal distribution.
        // (1-beta)*N(x, (2.38)^2 Cov / n_pars) + beta*N(x, (0.1)^2 I / n_pars)
        // Only start adapting the walk after a defined number of iterations.
        if (iter >= RVF_BEGIN_ADAPT)
        {
            // Calculate the mean and covariance of the chain every so often.
            if (iter % RVF_UPDATE_COV == 0)
            {
                ConstructMeanCov(iter);
            }
            
            // Construct the part of the proposal covariance based on the empirical covariance.
            gsl_matrix* emp_cov = gsl_matrix_alloc(n_pars, n_pars);
            for (int i = 0; i < n_pars; ++i)
            {
                for (int j = 0; j < n_pars; ++j)
                {
                    gsl_matrix_set(emp_cov, i, j, covariance[i][j]);
                }
            }
            
            // Scale the empirical covariance matrix appropriately.
            // Recall that scaling a multi-variate normal distribution by a constant, scales the
            // covariance by the square of the constant.
            gsl_matrix_scale(emp_cov, 2.38*2.38 / static_cast<double>(n_pars));
            
            // Scale the identity matrix (random walk) and add together.
            gsl_matrix_scale(cov, beta_prob*beta_prob);
            gsl_matrix_add(cov, emp_cov);
            gsl_matrix_free(emp_cov);
        }
        
        // Perform the Cholesky decomposition of the covariance matrix.
        gsl_linalg_cholesky_decomp(cov);

        // Define a vector of current values.
        gsl_vector* mu = gsl_vector_alloc(n_pars);
        for (int i = 0; i < n_pars; ++i)
        {
            gsl_vector_set(mu, i, 0.0);
        }

        // Define a result vector.
        gsl_vector* result = gsl_vector_alloc(n_pars);

        // Propose new parameter set.
        gsl_ran_multivariate_gaussian(rng, mu, cov, result);

        // Store the result in the next position in the chain.
        for (int i = 0; i < n_pars; ++i)
        {
            theta[i]->chain[iter + 1] = theta[i]->chain[iter] + gsl_vector_get(result, i);
        }
        
        // Clear the spaces used in calculating the covariance matrix of the proposal distribution.
        gsl_vector_free(result);
        gsl_vector_free(mu);
        gsl_matrix_free(cov);
    }
    else
    {
        // Randomly sample the ID of another MCMC chain.
        int alt_id;
        do
        {
            alt_id = gsl_rng_uniform_int(rng, RVF_CHAINS);
        } while (alt_id == id);

        // Declare space for the mean and covariance.
        gsl_vector* mu = gsl_vector_alloc(n_pars);
        gsl_matrix* cov = gsl_matrix_alloc(n_pars, n_pars);

        // Get that chains current mean and covariance.
        // Race conditions here... parameter set would be rejected, but
        // may get covariance matrix that is not invertible.
        for (int i = 0; i < n_pars; ++i)
        {
            gsl_vector_set(mu, i, mcmc[alt_id]->mean[i]);
            for (int j = 0; j < n_pars; ++j)
            {
                gsl_matrix_set(cov, i, j, mcmc[alt_id]->covariance[i][j]);
            }
        }

        // Define a result vector.
        gsl_vector* result = gsl_vector_alloc(n_pars);

        // Propose new parameter set.
        gsl_ran_multivariate_gaussian(rng, mu, cov, result);

        // Store the result in the next position in the chain.
        for (int i = 0; i < n_pars; ++i)
        {
            theta[i]->chain[iter + 1] = gsl_vector_get(result, i);
        }
        
        // Clear the spaces used in calculating the covariance matrix of the proposal distribution.
        gsl_vector_free(result);
        gsl_vector_free(mu);
        gsl_matrix_free(cov);
    }
}

// Check the validity of proposed parameters.
bool Mcmc::IsProposalValid(int iter)
{
    // Declare that all proposed parameters are initially valid.
    bool is_valid = true;

    // For each parameter, use member function to check it is valid.
    for (int i = 0; i < n_pars; ++i)
    {
        is_valid &= theta[i]->IsValid(iter);
    }

    // Return if all parameters are valid.
    return is_valid;
}

// Function to print the progress of the MCMC to console.
void Mcmc::PrintProgress(int iter)
{
    if (iter % RVF_UPDATE_FREQ == 0)
    {
        // Note how many iterations have been completed.
        std::cout << "\nIterations completed: " << iter << " / " << RVF_NMCMC << std::endl;

        // Note the acceptance rate.
        std::cout << "Acceptance rate: " << round(acceptance_rate * 100 * 10) / 10.0;
        std::cout << "%" << std::endl;

        // Note the current log-likelihood.
        std::cout << "Log-likelihood: " << log_likelihood_chain[iter] << std::endl;

        // Note the current parameter values.
        for (int i = 0; i < n_pars; ++i)
        {
            std::cout << theta[i]->par_name << ": " << theta[i]->chain[iter] << std::endl;
        }        
    }
}

// Function to evaluate the prior at a given iteration in the chain.
double Mcmc::LogPrior(int iter)
{
    // Initialise the log-prior.
    double log_prior = 0.0;

    // Assuming all parameters are independent, for each parameter...
    for (int i = 0; i < n_pars; ++i)
    {
        // Evaluate the prior at the given iteration.
        log_prior += theta[i]->Lpdf(iter);
    }

    // Return the log-prior ratio.
    return log_prior;
}

// Function to write priors to file.
void Mcmc::WritePriors()
{
    // Delcare the file stream.
    std::fstream file;

    // Open the file for writing.
    file.open(static_cast<std::string>(RVF_ODIR) + "mcmc_priors.csv", std::fstream::out);

    // Check that the file opened.
    if (file.is_open())
    {
        // If the file could be opened, write the header.
        file << "PAR_NAME,PRIOR_DIST,PRIOR_PARS";
        
        // For each parameter that we are inferring, output its information.
        for (int i = 0; i < n_pars; ++i)
        {
            file << "\n" << theta[i]->GetInfo();
        }

        // Close the file once completed.
        file.close();
    }
    else
    {
        std::cout << "Prior file could not be opened for writing." << std::endl;
        exit(1);
    }    
}

// Function to write the chains of the inferred parameters to file.
void Mcmc::WriteChains(int iter)
{
    // Declare the file stream.
    std::fstream file;

    // Define the start and end iteration to save the chains for.
    int start_iter = ((iter-1) / RVF_UPDATE_FREQ)*RVF_UPDATE_FREQ;
    int end_iter = iter;

    // If we're at the start, open the file without appending.
    if (start_iter == 0)
    {
        // Keep trying to open the file until successful. (infinite loop warning :))
        do {
            file.open(static_cast<std::string>(RVF_ODIR) + "mcmc_" + std::to_string(id) + ".csv", std::fstream::out);
        } while (!file.is_open());

        // Also output headers to the file!
        file << "STEP,PAR_NAME,PAR_VALUE,LOG_LIKELIHOOD,LOG_PRIOR";
    }
    else
    {
        // Otherwise, open with appending.
        // Keep trying to open the file until successful. (infinite loop warning :))
        do {
            file.open(static_cast<std::string>(RVF_ODIR) + "mcmc_" + std::to_string(id) + ".csv", std::fstream::app);
        } while (!file.is_open());
    }

    // For each step in the MCMC iteration, output
    // the name of the parameter and the value in the chain.
    for (int step = start_iter; step < end_iter; ++step)
    {
        for (int i = 0; i < n_pars; ++i)
        {
            file << "\n" << std::fixed << std::setprecision(3);
            file << step << "," << theta[i]->par_name;
            file << "," << theta[i]->chain[step];
            file << "," << log_likelihood_chain[step];
            file << "," << log_prior_chain[step];
        }
    }

    // Close the file after writing has finished.
    file.close();    
};
