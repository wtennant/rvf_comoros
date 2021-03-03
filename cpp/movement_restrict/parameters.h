 // parameters.h: Defines the class named Parameters which
// are used with in the model.
#include <string>   // Using strings.

#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters
{
public:
    int n_age;                  // Number of livestock age-groups (grouped by year).
    int n_age_move;             // Maximum age group which can be moved.
    int n_island;               // Nimber of islands.
    int n_steps;                // Number of time steps to simulate (epi-weeks).
    double import_start;        // Start time of import into Grande Comore.
    double import_duration;     // Duration of import into Grande Comore.
    double import_size;         // Number of infectious imports into Grande Comore per year.
    double import_start_2;      // Start time of import into Grande Comore.
    double import_duration_2;   // Duration of import into Grande Comore.
    double import_size_2;       // Number of infectious imports into Grande Comore per year.
    double* n_pop;              // Population size of each island.
    double* p_immune_start;     // Proportion of animals that are immune at time zero on each island.
    double* trans_scale;        // Scalar for transmission rate on each island.
    double* ndvi_rate;          // Scalar for how NDVI changes force of infection for each island.
    double** move;              // Number of annual movements between each island.
    double* pop_structure;      // Population structure (proportion of population in
                                // each age group).
    double* mortality;          // Age-specific mortality rates (proportion per epi-week).

    // Constructor and destructor.
    Parameters();
    ~Parameters();

private:
    void WriteDefaultPars();    // Write default parameters to file.
};

#endif