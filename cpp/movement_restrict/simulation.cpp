// simulation.cpp: Definition of the function that simulates RVF.
#include <chrono>               // Seed for random number generation.
#include <fstream>              // Writing data to file.
#include <iostream>             // Input-output to console.
#include <iomanip>
#include <math.h>               // Mathematical operations.
#include "omp.h"                // Parallelism.
#include "config.h"             // Location to output file.
#include "Simulation.h"         // Definition of the Simulation class.

// Constructor to allocate memory for each infection state, age and time step.
Simulation::Simulation(Parameters* pars)
{
    // Store local copies of the number of age groups, islands and time steps.
    n_age = pars->n_age;
    n_age_move = pars->n_age_move;
    n_island = pars->n_island;
    n_steps = pars->n_steps;

    // Calculate the scaling factor for imports in each age group.
    imp_scaling = 0.0;
    for (int age = 0; age < n_age_move; ++age)
    {
        imp_scaling += pars->pop_structure[age]; 
    }    

    // Allocate space for each time step. 
    S = new double**[n_steps + 1];
    E = new double**[n_steps + 1];
    I = new double**[n_steps + 1];
    R = new double**[n_steps + 1];
    exports = new double**[n_steps + 1];
    imports = new double**[n_steps + 1];

    // For each time step, allocate space for each island.
    for (int t = 0; t <= n_steps; ++t)
    {
        S[t] = new double*[n_island];
        E[t] = new double*[n_island];
        I[t] = new double*[n_island];
        R[t] = new double*[n_island];
        exports[t] = new double*[n_island];
        imports[t] = new double*[n_island];
            
        // For each island, allocate age groups.
        for (int i = 0; i < n_island; ++i)
        {
            S[t][i] = new double[n_age];
            E[t][i] = new double[n_age];
            I[t][i] = new double[n_age];
            R[t][i] = new double[n_age];
            exports[t][i] = new double[n_age];
            imports[t][i] = new double[n_age];

            for (int age = 0; age < n_age; ++age)
            {
                S[t][i][age] = 0.0;
                E[t][i][age] = 0.0;
                I[t][i][age] = 0.0;
                R[t][i][age] = 0.0;
                exports[t][i][age] = 0.0;
                imports[t][i][age] = 0.0;
            }
        }
    }

    // Allocate space for the number of movements of animals between islands at
    // each time step.
    move_S = new double**[n_island];
    move_E = new double**[n_island];
    move_I = new double**[n_island];
    move_R = new double**[n_island];
    for (int i = 0; i < n_island; ++i)
    {
        move_S[i] = new double*[n_island];
        move_E[i] = new double*[n_island];
        move_I[i] = new double*[n_island];
        move_R[i] = new double*[n_island];
        for (int j = 0; j < n_island; ++j)
        {
            move_S[i][j] = new double[n_age_move];
            move_E[i][j] = new double[n_age_move];
            move_I[i][j] = new double[n_age_move];
            move_R[i][j] = new double[n_age_move];
            for (int age = 0; age < n_age_move; ++age)
            {
                move_S[i][j][age] = 0.0;
                move_E[i][j][age] = 0.0;
                move_I[i][j][age] = 0.0;
                move_R[i][j][age] = 0.0;
            }
        }
    }

    // Open up a file stream for recording simulated data.
    file.open(static_cast<std::string>(RVF_ODIR) + "simulation.csv", std::fstream::out);

    // Check that the file is ready for writing.
    if (!file.is_open())
    {
        std::cout << "Error: there was a problem opening the file for writing";
        std::cout << " simulation output." << std::endl;
        exit(1);
    }
    else
    {
        // Write headers to file.
        file << "SAMPLE_ID,TIME,ISLAND_ID,AGE,S,E,I,R,IMPORTS_EI,EXPORTS_EI";
    }
}

// Destructor to de-allocate memory for each infection state, age and time step.
Simulation::~Simulation()
{
    // Free-up memory that was allocated for movement.
    for (int i = 0; i < n_island; ++i)
    {
        for (int j = 0; j < n_island; ++j)
        {
            delete[] move_R[i][j];
            delete[] move_I[i][j];
            delete[] move_E[i][j];
            delete[] move_S[i][j];
        }
        delete[] move_R[i];
        delete[] move_E[i];
        delete[] move_I[i];
        delete[] move_S[i];
    }
    delete[] move_R;
    delete[] move_I;
    delete[] move_E;
    delete[] move_S;

    // Free-up memory that was allocated for the compartments.
    for (int t = 0; t <= n_steps; ++t)
    {
        for (int i = 0;i < n_island; ++i)
        {
            // Free-up memory that was allocated for each age.
            delete[] S[t][i];
            delete[] E[t][i];
            delete[] I[t][i];
            delete[] R[t][i];
            delete[] exports[t][i];
            delete[] imports[t][i];
        }

        // Free up memory that was allocated for each island.
        delete[] S[t];
        delete[] E[t];
        delete[] I[t];
        delete[] R[t];
        delete[] exports[t];
        delete[] imports[t];
    }

    // Free-up the memory that was allocated for each time step.
    delete[] S;
    delete[] E;
    delete[] I;
    delete[] R;
    delete[] exports;
    delete[] imports;

    // Close the file used to record simulation output.
    file.close();
}

// Function to simulate the RVF given model parameters and NDVI data.
void Simulation::Simulate(Data* data, Parameters* pars, int start_time, int end_time)
{
    // If at the first time point, initialise the compartments.
    if (start_time == 0)
    {
        // For each island in the metapopulation...
        for (int i = 0; i < n_island; ++i)
        {
            // Sample the number of individuals that are immune at the start of the simulation.
            double init_immune = pars->p_immune_start[i]*pars->n_pop[i];

            // Sample the age distriubtion of the initial livestock populations from a multinomial distribution.
            for (int age = 0; age < n_age; ++age)
            {
                S[0][i][age] = pars->pop_structure[age]*(pars->n_pop[i] - 10.0 - init_immune);
                E[0][i][age] = pars->pop_structure[age]*5.0;
                I[0][i][age] = pars->pop_structure[age]*5.0;
                R[0][i][age] = pars->pop_structure[age]*init_immune;
            }
        }
    }

    // Calculate the probability of the incubation and infectious period ending.
    // This saves the calculation of the exponential function an excessive
    // number of times.
    double prob_inc_end = 1.0;
    double prob_inf_end = 1.0;

    // Create a counter to keep track of the number of individuals on each island
    // and the number of infectious individuals on each island.
    double* N = new double[n_island]{0.0};
    double* total_infectious = new double[n_island]{0.0};
    double* deaths = new double[n_island]{0.0};

    // Keep track of the net movement numbers on each island.
    double* net_move = new double[n_island]{0.0};

    // Run the model for the desired number of time steps.
    for (int t = start_time; t < end_time; ++t)
    {
        // For each island, run the infection and ageing process.
        for (int i = 0; i < n_island; ++i)
        {
            // Reset the total number of individuals counters.
            N[i] = 0.0;
            total_infectious[i] = 0.0; 
            deaths[i] = 0.0;
            net_move[i] = 0.0;

            // Calculate the number of deaths in each age group.
            for (int age = 0; age < n_age; ++age)
            {
                S[t + 1][i][age] = (1.0 - pars->mortality[age])*S[t][i][age];
                E[t + 1][i][age] = (1.0 - pars->mortality[age])*E[t][i][age];
                I[t + 1][i][age] = (1.0 - pars->mortality[age])*I[t][i][age];
                R[t + 1][i][age] = (1.0 - pars->mortality[age])*R[t][i][age];

                // Count up the total number of individuals and infectious individuals.
                total_infectious[i] += I[t][i][age];
                N[i] += S[t][i][age] + E[t][i][age] + I[t][i][age] + R[t][i][age];
                deaths[i] += pars->mortality[age]*(S[t][i][age] + E[t][i][age] + I[t][i][age] + R[t][i][age]);

                // Reset the total number of imports and exports.
                imports[t + 1][i][age] = 0.0;
                exports[t + 1][i][age] = 0.0;
            }

            // Define the NDVI minimum to be used when calculating the transmission rate.
            double min_ndvi;
            switch(RVF_MIN_NDVI)
            {
                case 0:     // Use minimum of each island.
                    min_ndvi = data->ndvi.min_local_ndvi[i];
                    break;
                case 1:     // Use the minimum across all islands.
                    min_ndvi = data->ndvi.min_global_ndvi;
                    break;
                case 2:     // Don't use the minimum.
                    min_ndvi = 0.0;
                    break;
                default:    // By default, use the global minimum.
                    min_ndvi = data->ndvi.min_global_ndvi;
            }

            // Calculate the force of infection.
            double foi;

            // First calculate the influence of ndvi on the transmission rate.
            // Define which model to use.
            switch(RVF_TRANSMISSION)
            {
                case 0:     // Constant transmission.
                    foi = pars->trans_scale[i];
                case 1:     // Linear transmission model.
                    foi = pars->ndvi_rate[i]*(data->ndvi.ndvi[i][t] - min_ndvi) + pars->trans_scale[i];
                case 2:     // Exponential transmission model.
                    foi = exp(pars->ndvi_rate[i]*(data->ndvi.ndvi[i][t] - min_ndvi) + pars->trans_scale[i]);
            }

            // Then calculate the probability of becoming infected assuming a Poisson process.
            foi = 1.0 - exp(-foi*total_infectious[i] / N[i]);
            
            // Transition individuals between infection classes.
            for (int age = 0; age < n_age; ++age)
            {
                // Sample the number of individuals finishing their incubation
                // and infectious periods.
                double inf_end = prob_inf_end*I[t + 1][i][age];
                double inc_end = prob_inc_end*E[t + 1][i][age];

                // Sample the number of infections.
                double infs = foi*S[t + 1][i][age];
                
                // Carry out the transitions.
                // As we're working with double, safer to do additions first.
                R[t + 1][i][age] += inf_end;
                I[t + 1][i][age] += inc_end;
                I[t + 1][i][age] -= inf_end;
                E[t + 1][i][age] += infs;
                E[t + 1][i][age] -= inc_end;
                S[t + 1][i][age] -= infs;

                // Update the total number of infectious individuals on the island.
                total_infectious[i] += inc_end;
                total_infectious[i] -= inf_end;
            }

            // Transition individuals between age groups.
            for (int age = 0; age < n_age - 1; ++age)
            {
                // Sample the number of individuals which age.
                double S_age = S[t + 1][i][age] / 12.0 / 4.0;
                double E_age = E[t + 1][i][age] / 12.0 / 4.0;
                double I_age = I[t + 1][i][age] / 12.0 / 4.0;
                double R_age = R[t + 1][i][age] / 12.0 / 4.0;

                // Add individuals who survived into the next age group.
                S[t + 1][i][age + 1] += S_age;
                E[t + 1][i][age + 1] += E_age;
                I[t + 1][i][age + 1] += I_age;
                R[t + 1][i][age + 1] += R_age;

                // Take off the individuals who survived from the current age group.
                S[t + 1][i][age] -= S_age;
                E[t + 1][i][age] -= E_age;
                I[t + 1][i][age] -= I_age;
                R[t + 1][i][age] -= R_age;
            }
        }

        // Calculate the number of movements between islands.
        CalculateMove(data, pars, t+1);

        // Carry out the movements between each island for all age groups.
        for (int i = 0; i < n_island; ++i)
        {
            for (int j = 0; j < n_island; ++j)
            {
                for (int age = 0; age < n_age_move; ++age)
                {
                    // Move all animals between the islands.
                    // Do this for each compartment in the infection process.
                    S[t + 1][j][age] += move_S[i][j][age];
                    S[t + 1][i][age] -= move_S[i][j][age];
                    E[t + 1][j][age] += move_E[i][j][age];
                    E[t + 1][i][age] -= move_E[i][j][age];
                    I[t + 1][j][age] += move_I[i][j][age];
                    I[t + 1][i][age] -= move_I[i][j][age];
                    R[t + 1][j][age] += move_R[i][j][age];
                    R[t + 1][i][age] -= move_R[i][j][age];

                    // Calculate the net movements between island i and j.
                    net_move[i] -= (move_S[i][j][age] + move_E[i][j][age] +
                        move_I[i][j][age] + move_R[i][j][age]);
                    net_move[j] += (move_S[i][j][age] + move_E[i][j][age] +
                        move_I[i][j][age] + move_R[i][j][age]);

                    // Get all the imports and exports between the islands.
                    exports[t + 1][i][age] += (move_E[i][j][age] + move_I[i][j][age]);
                    imports[t + 1][j][age] += (move_E[i][j][age] + move_I[i][j][age]);
                }
            }
        }
        
        // Carry out movement imports if in the appropriate time window.
        if ((t >= pars->import_start) && (t < pars->import_start + pars->import_duration))
        {               
            // Import infectious animals into Grande Comore.
            for (int age = 0; age < n_age_move; ++age)
            {
                I[t + 1][1][age] += pars->import_size*pars->pop_structure[age] / imp_scaling;
                net_move[1] += pars->import_size*pars->pop_structure[age] / imp_scaling;
                //imports[t + 1][1][age] += pars->import_size*pars->pop_structure[age] / imp_scaling;
            }
        }

        // Carry out movement imports from a second introduction.
        switch (RVF_SEC_IMPORT)
        {
            case 0:
                break;
            case 1:
                if ((t >= pars->import_start_2) && (t < pars->import_start_2 + pars->import_duration_2))
                {               
                    // Import infectious animals into Grande Comore.
                    for (int age = 0; age < n_age_move; ++age)
                    {
                        I[t + 1][1][age] += pars->import_size_2*pars->pop_structure[age] / imp_scaling;
                        net_move[1] += pars->import_size_2*pars->pop_structure[age] / imp_scaling;
                        //imports[t + 1][1][age] += pars->import_size_2*pars->pop_structure[age] / imp_scaling;
                    }
                }
                break;
            default:
                break;
        }

        // For every island...
        for (int i = 0; i < n_island; ++i)
        {
            // ...rebalance the total number of individuals on island i.
            S[t + 1][i][0] += (deaths[i] - net_move[i]);
        }
    }
                  
    // De-allocate memory for the total count and infectious count.
    delete[] N;
    delete[] total_infectious;
}

// Function to calculate the number of movements between islands at
// a given time step.
void Simulation::CalculateMove(Data* data,
                               Parameters* pars,
                               const int t)
{
    // Set up counter to store the number of individuals that are up for
    // moving on each island.
    double* total_move_N = new double[n_island]{0.0};

    // Calculate total number of individuals who can be moved on each island.
    for (int i = 0; i < n_island; ++i)
    {
        for (int age = 0; age < n_age_move; ++age)
        {
            total_move_N[i] += S[t][i][age] + E[t][i][age] +
                I[t][i][age] + R[t][i][age];
        }
    }

    // Move individuals from island i...
    for (int i = 0; i < n_island; ++i)
    {
        // Calculate the number of animals to pass to other islands.
        for (int j = 0; j < n_island; ++j)
        {
            // If there are movements between island i and j, then...
            if ((i != j) && (pars->move[i][j] > 0.0))
            {
                // ...calculate the total number of movements between island i and island j.
                double movements = pars->move[i][j];

                // For each age group, transfer the infected individuals from island i to island j.
                for (int age = 0; age < n_age_move; ++age)
                {
                    move_S[i][j][age] = movements * (S[t][i][age] / total_move_N[i]);
                    move_E[i][j][age] = movements * (E[t][i][age] / total_move_N[i]);
                    move_I[i][j][age] = movements * (I[t][i][age] / total_move_N[i]);
                    move_R[i][j][age] = movements * (R[t][i][age] / total_move_N[i]);
                }
            }
        }
    }

    // Free up memory for the total number of individuals which can be moved.
    delete[] total_move_N;
}

// Function to write simulation output to file.
void Simulation::WriteOutput(int iter)
{
    // For every time step...
    for (int t = 0; t <= n_steps; ++t)
    {
        // For every island...
        for (int i = 0; i < n_island; ++i)
        {
            // ...for every age group...
            for (int age = 0; age < n_age; ++age)
            {
                // Write the number of individuals in each infection compartment
                // to the file.
                file << "\n" << std::fixed << std::setprecision(3);
                file << iter << "," << t << "," << i << "," << age << "," << S[t][i][age];
                file << "," << E[t][i][age] << "," << I[t][i][age] << "," << R[t][i][age];
                file << "," << imports[t][i][age] << "," << exports[t][i][age];
            }
        }
    }

    // Send the contents of the filestream to the file.
    file.flush();
}