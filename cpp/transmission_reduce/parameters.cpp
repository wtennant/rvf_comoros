// parameters.cpp: Contains the constructor and destructor of the
// class, Parameters.
#include <fstream>          // File stream.
#include <iostream>         // Input-output stream.
#include "config.h"         // Output directory.
#include "parameters.h"     // Definition of the Parameters class.

// Default constructor for the Parameters class.
Parameters::Parameters()
{
    // Number of livestock age-groups (grouped by year).
    // Note: Currently undefined behaviour if n_age is not equal to 10.
    n_age = 10;

    // Maximum age group which can be moved.
    n_age_move = 2;

    // Define the number of islands in the metapopulation model.
    n_island = 4;

    // Number of livestock in the population.
    n_pop = new double[n_island];
    n_pop[0] = 93616.0;     // Anjouan.
    n_pop[1] = 224353.0;    // Grande Comore.
    n_pop[2] = 20052.0;     // Mayotte.
    n_pop[3] = 31872.0;     // Moheli.

    // Number of time steps to simulate (epi-weeks).
    n_steps = 576;

    // Proportion of animals that are immune at time zero on each island.
    p_immune_start = new double[n_island]{0.1298};

    // Scalars for how NDVI influences the force of infection (FoI).
    // FoI = 1 - exp(-exp(-rate*NDVI + constant)*I/N).
    trans_scale = new double[n_island]{-2.19};
    ndvi_rate = new double[n_island]{3.15};

    // Vector control level on each island.
    vec_control = new double[n_island]{0.0};
    for (int i = 0; i < n_island; ++i)
    {
        vec_control[i] = 0.0;
    }

    // Declare space for the number of movements between islands.
    move = new double*[n_island];
    for (int i = 0; i < n_island; ++i)
    {
        move[i] = new double[n_island]{0.0};
    }

    // Declare space for each of the age-dependent parameters.
    // Proportion of the population in each age group.
    pop_structure = new double[n_age]{0.294, 0.206, 0.144, 0.101, 0.071,
                                      0.050, 0.035, 0.024, 0.017, 0.058};

    // Mortality rate for each age group.
    mortality = new double[n_age]{0.0088, 0.0088, 0.0088, 0.0088, 0.0088,
                                  0.0088, 0.0088, 0.0088, 0.0088, 0.0062};

    // Write the parameters to file. Only do this once.
    WriteDefaultPars();
}

// Default deconstructor for the Parameter class.
Parameters::~Parameters()
{
    // Free dynamically allocated memory.
    for (int i = 0; i < n_island; ++i)
    {
        delete[] move[i];
    }
    delete[] move;
    delete[] vec_control;
    delete[] trans_scale;
    delete[] p_immune_start;
    delete[] mortality;
    delete[] pop_structure;
    delete[] n_pop;
}

// Write default parameters to file.
void Parameters::WriteDefaultPars()
{
    // Declare the file stream.
    std::fstream file;

    // Open the file to record parameters.
    file.open(static_cast<std::string>(RVF_ODIR) + "default_pars.csv", std::fstream::out);

    // Check that the file is open.
    if (file.is_open())
    {
        // Record the header for the default parameter file.
        file << "PAR_NAME,PAR_VALUE" << std::endl;

        // For each parameter, output the name and value.
        file << "n_age," << n_age << std::endl;
        file << "nsteps," << n_steps << std::endl;
        file << "n_pop,[";
        for (int i = 0; i < n_island; ++i)
        {
            file << n_pop[i] << ((i < n_island - 1) ? " " : "");
        }
        file << "]" << std::endl;
        file << "pop_structure,[";
        for (int age = 0; age < n_age; ++age)
        {
            file << pop_structure[age] << ((age < n_age - 1) ? " " : "");
        }
        file << "]" << std::endl;
        file << "mortality,[";
        for (int age = 0; age < n_age; ++age)
        {
            file << mortality[age] << ((age < n_age - 1) ? " " : "");
        }
        file << "]" << std::endl;
    }
    else
    {
        // If the file could not be opened, throw a wobbly.
        std::cout << "File could not be opened to record default parameter set." << std::endl;
        exit(1);
    }

    // Close the file after writing has finished.
    file.close();    
}