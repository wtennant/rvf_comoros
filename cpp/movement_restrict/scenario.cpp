// scenario.cpp: Defines the function to run all the scenario tests.
#include <chrono>               // Time used as the random number generator seed.
#include <fstream>              // Writing to file.
#include <iostream>             // Output information to console.
#include <iomanip>              // Set precision of file output.
#include "gsl/gsl_randist.h"    // Distributions to evaluate likelihood.
#include "config.h"             // Configuration of file output.
#include "scenario.h"           // Definition of the Scenario class.

// Constructor of the Scenario class.
Scenario::Scenario()
{
    // Define the number of parameters involed in the scenario testing.
    n_scenarios = 0;    // The number of parameters independent of scenario choice.

    // Add in the number of parameters depending on the scenarios turned on.
    switch(RVF_SCALE_MAY)
    {
        case 0:     // Movement into Mayotte directly from posterior samples.
            break;
        case 1:     // Movement into Mayotte is scaled by scale_move_anj_mayotte.
            ++n_scenarios;
            break;
        default:    // Default is always the scenario off.
            break;
    }

    // Scenario with a second importation window into Grande Comore.
    switch(RVF_SEC_IMPORT)
    {
        case 0:
            break;
        case 1:
            n_scenarios = n_scenarios + 3;
            break;
        default:
            break;
    }

    // Set up the random number generator based on time.
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, 2*std::chrono::system_clock::now().time_since_epoch().count());

    // Space to store all the scenario parameters.
    theta = new Theta*[n_scenarios];

    // Define an index which keeps track of how many scenarios have been stored.
    int scenario_idx = 0;

    // What scenario parameters are in there depends upon the scenarios turned on.
    switch(RVF_SCALE_MAY)
    {
        case 0:     // Movement into Mayotte directly from posterior samples.
            break;
        case 1:     // Movement into Mayotte is scaled by scale_move_anj_mayotte.
            theta[scenario_idx] = &scale_move_anj_may;
            ++scenario_idx;
            break;
        default:    // Default is always the scenario off.
            break;
    }
    switch(RVF_SEC_IMPORT)
    {
        case 0:     // No second importation into Grande Comore.
            break;
        case 1:     // Second importation into Grande Comore
            theta[scenario_idx] = &import_duration_two;
            theta[scenario_idx + 1] = &import_size_two;
            theta[scenario_idx + 2] = &import_start_two;
            scenario_idx = scenario_idx + 3;
            break;
        default:
            break;
    }
    
    // For each scenario parameter, allocate space for its samples.
    for (int i = 0; i < n_scenarios; ++i)
    {
        theta[i]->chain = new double[RVF_NSAMPLES];
    }

    // Declare space for the IDs of posterior samples.
    post_idx = new int[RVF_NSAMPLES];
}

// Deconstructor of the Scenario class.
Scenario::~Scenario()
{
    // De-allocate the random number generator.
    gsl_rng_free(rng);

    // Free up the scenario parameters.
    for (int i = 0; i < n_scenarios; ++i)
    {
        delete[] theta[i]->chain;
    }
    delete[] theta;

    // Free up the posterior sample IDs.
    delete[] post_idx;
}

// Scenario testing algorithm:
// Randomly sample from the posterior distribution, then randomly sample 
// from the distribution of each scenario parameter, then run the simulation,
// then store the result.
void Scenario::Run(Data* data, Parameters* pars, Posterior* post, Simulation* sim)
{
    // For each requested scenario sample:
    for (int i = 0; i < RVF_NSAMPLES; ++i)
    {
        // Randomly sample which posterior sample to use.
        int post_sample_id = gsl_rng_uniform_int(rng, post->n_samples);
        post_idx[i] = post_sample_id;

        // Set the simulation parameters from this posterior sample.
        SetSimPars(post_sample_id, pars, post);

        // Randomly sample and store all scenario parameters.
        for (int j = 0; j < n_scenarios; ++j)
        {
            theta[j]->Sample(i);
        }

        // Update the simulation parameters based on sampled scenario parameters.
        UpdateSimPars(i, pars);

        // Run the simulation for all time steps.
        sim->Simulate(data, pars, 0, pars->n_steps);

        // Save the simulation result.
        sim->WriteOutput(i);

        // Update the user on progress of the scenario tests.
        PrintProgress(i);
    }

    // Do the final update of scenario tests completed.
    PrintProgress(RVF_NSAMPLES);

    // Write the sampled scenario parameters to file.
    WriteScenarioPars();

    // Write the sampled posterior values to file.
    WritePosteriorSamples(post);
}

// Function to set the simulation parameters from the posterior.
void Scenario::SetSimPars(int sample_id, Parameters* pars, Posterior* post)
{
    // Set the parameters in the simulation.
    pars->p_immune_start[0] = post->p_immune_start_anj.chain[sample_id];
    pars->p_immune_start[1] = post->p_immune_start_gra.chain[sample_id];
    pars->p_immune_start[2] = post->p_immune_start_may.chain[sample_id];
    pars->p_immune_start[3] = post->p_immune_start_moh.chain[sample_id];
    pars->move[0][1] = post->move_anj_gra.chain[sample_id];
    pars->move[0][2] = post->move_anj_may.chain[sample_id];
    pars->move[0][3] = post->move_anj_moh.chain[sample_id];
    pars->move[1][0] = post->move_gra_anj.chain[sample_id];
    pars->move[1][3] = post->move_gra_moh.chain[sample_id];
    pars->move[3][0] = post->move_moh_anj.chain[sample_id];
    pars->move[3][1] = post->move_moh_gra.chain[sample_id];
    pars->import_start = post->import_start.chain[sample_id];
    pars->import_duration = post->import_duration.chain[sample_id];
    pars->import_size = post->import_size.chain[sample_id];

    // The model parameters used depend on configuration.
    switch(RVF_DIFF_NDVI)
    {
        case 0:     // Same NDVI per island.
            pars->ndvi_rate[0] = post->ndvi_rate.chain[sample_id];
            pars->ndvi_rate[1] = post->ndvi_rate.chain[sample_id];
            pars->ndvi_rate[2] = post->ndvi_rate.chain[sample_id];
            pars->ndvi_rate[3] = post->ndvi_rate.chain[sample_id];
            break;
        case 1:     // Different NDVI per island.
            pars->ndvi_rate[0] = post->ndvi_rate_anj.chain[sample_id];
            pars->ndvi_rate[1] = post->ndvi_rate_gra.chain[sample_id];
            pars->ndvi_rate[2] = post->ndvi_rate_may.chain[sample_id];
            pars->ndvi_rate[3] = post->ndvi_rate_moh.chain[sample_id];
            break;
        default:    // Default is same.
            pars->ndvi_rate[0] = post->ndvi_rate.chain[sample_id];
            pars->ndvi_rate[1] = post->ndvi_rate.chain[sample_id];
            pars->ndvi_rate[2] = post->ndvi_rate.chain[sample_id];
            pars->ndvi_rate[3] = post->ndvi_rate.chain[sample_id];
    }
    switch(RVF_DIFF_SCALE)
    {
        case 0:     // Same transmission scale per island.
            pars->trans_scale[0] = post->trans_scale.chain[sample_id];
            pars->trans_scale[1] = post->trans_scale.chain[sample_id];
            pars->trans_scale[2] = post->trans_scale.chain[sample_id];
            pars->trans_scale[3] = post->trans_scale.chain[sample_id];
            break;
        case 1:     // Different transmission scale per island.
            pars->trans_scale[0] = post->trans_scale_anj.chain[sample_id];
            pars->trans_scale[1] = post->trans_scale_gra.chain[sample_id];
            pars->trans_scale[2] = post->trans_scale_may.chain[sample_id];
            pars->trans_scale[3] = post->trans_scale_moh.chain[sample_id];
            break;
        default:    // Default is same.
            pars->trans_scale[0] = post->trans_scale.chain[sample_id];
            pars->trans_scale[1] = post->trans_scale.chain[sample_id];
            pars->trans_scale[2] = post->trans_scale.chain[sample_id];
            pars->trans_scale[3] = post->trans_scale.chain[sample_id];
    }
}

// Function to update the simulation parameters based on scenario tests.
void Scenario::UpdateSimPars(int sample_id, Parameters* pars)
{
    // Parameters to update depends on the secnarios requested.
    switch(RVF_SCALE_MAY)
    {
        case 0:     // Movement into Mayotte directly from posterior samples.
            break;
        case 1:     // Movement into Mayotte is scaled by scale_move_anj_may.
            pars->move[0][2] *= scale_move_anj_may.chain[sample_id];
            break;
        default:    // Default is always the scenario off.
            break;
    }
    
    // Second importation into Grande Comore.
    switch(RVF_SEC_IMPORT)
    {
        case 0:     // No second importation.
            break;
        case 1:     // Allow second importation into Grande Comore.
            pars->import_duration_2 = import_duration_two.chain[sample_id];
            pars->import_size_2 = import_size_two.chain[sample_id];
            pars->import_start_2 = import_start_two.chain[sample_id];
        default:    // No second importation.
            break;
    }
    
    // Isolate Grande Comore from the network.
    switch(RVF_ISOLATE_GRA)
    {
        case 0:     // Use existing network structure.
            break;
        case 1:     // Set all movements into and out of Grande Comore to zero.
            pars->move[1][0] *= 0.0;
            pars->move[1][3] *= 0.0;
            pars->move[0][1] *= 0.0;
            pars->move[3][1] *= 0.0;
            break;
        default:    // Default to existing network structure.
            break;
    }

    // Isolate Moheli from the network.
    switch(RVF_ISOLATE_MOH)
    {
        case 0:     // Use existing network structure.
            break;
        case 1:     // Set all movements into and out of Moheli to zero.
            pars->move[3][0] *= 0.0;
            pars->move[1][3] *= 0.0;
            pars->move[0][3] *= 0.0;
            pars->move[3][1] *= 0.0;
            break;
        default:    // Default to existing network structure.
            break;
    }
    
    // Isolate Anjouan from the network.
    switch(RVF_ISOLATE_ANJ)
    {
        case 0:     // Use existing network structure.
            break;
        case 1:     // Set all movements into and out of Anjouan to zero.
            pars->move[3][0] *= 0.0;
            pars->move[1][0] *= 0.0;
            pars->move[0][3] *= 0.0;
            pars->move[0][2] *= 0.0;
            pars->move[0][1] *= 0.0;
            break;
        default:    // Default to existing network structure.
            break;
    }
    
    // Isolate Mayotte from the network.
    switch(RVF_ISOLATE_MAY)
    {
        case 0:     // Use existing network structure.
            break;
        case 1:     // Set all movements into and out of Mayotte to zero.
            pars->move[0][2] *= 0.0;
            break;
        default:    // Default to existing network structure.
            break;
    }
}

// Function to print the progress of the scenario testing to console.
void Scenario::PrintProgress(int i)
{
    if (i % RVF_UPDATE_FREQ == 0)
    {
        // Note how many iterations have been completed.
        std::cout << "\rScenarios completed: " << i << " / " << RVF_NSAMPLES;
    }
}

// Function to write the sampled posterior parameters to file.
void Scenario::WritePosteriorSamples(Posterior* post)
{
    // Declare the file stream.
    std::fstream file;

    // Open a file for writing.
    file.open(static_cast<std::string>(RVF_ODIR) + "posterior.csv", std::fstream::out);

    // If we're at the start, open the file without appending.
    if (file.is_open())
    {
        // Output headers to the file!
        file << "SAMPLE_ID,PAR_NAME,PAR_VALUE";

        // For each scenario tested, save its ID, the scenario parameter names
        // and their values.
        for (int i = 0; i < RVF_NSAMPLES; ++i)
        {
            for (int j = 0; j < post->n_pars; ++j)
            {
                file << "\n" << std::fixed << std::setprecision(3);
                file << i << "," << post->theta[j]->par_name;
                file << "," << post->theta[j]->chain[post_idx[i]];
            }
        }
    }
    else
    {
        std::cout << "Error: there was a problem opening the file for writing";
        std::cout << " scenario output." << std::endl;
        exit(1);
    }

    // Close the file after writing has finished.
    file.close();    
};

// Function to write the scenario parameters to file.
void Scenario::WriteScenarioPars()
{
    // Declare the file stream.
    std::fstream file;

    // Open a file for writing.
    file.open(static_cast<std::string>(RVF_ODIR) + "scenario_pars.csv", std::fstream::out);

    // If we're at the start, open the file without appending.
    if (file.is_open())
    {
        // Output headers to the file!
        file << "SAMPLE_ID,PAR_NAME,PAR_VALUE";

        // For each scenario tested, save its ID, the scenario parameter names
        // and their values.
        for (int i = 0; i < RVF_NSAMPLES; ++i)
        {
            for (int j = 0; j < n_scenarios; ++j)
            {
                file << "\n" << std::fixed << std::setprecision(3);
                file << i << "," << theta[j]->par_name;
                file << "," << theta[j]->chain[i];
            }
        }
    }
    else
    {
        std::cout << "Error: there was a problem opening the file for writing";
        std::cout << " scenario output." << std::endl;
        exit(1);
    }

    // Close the file after writing has finished.
    file.close();    
};
