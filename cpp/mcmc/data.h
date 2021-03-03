// data.h: Defines the Seroprev and Data classes. These class contains the data
// needed to run and fit the RVF simulation model.
#include "parameters.h" // Definition of the Parameter class.

#ifndef DATA_H
#define DATA_H

class Ndvi
{
public:
    int n_island;                       // The number of islands.
    double start_month;                 // Starting month of the NDVI data.
    double start_year;                  // Starting year of the NDVI data.
    double** ndvi;                      // Time-series of the NDVI data.
    double* min_local_ndvi;             // Get the minimum NDVI on each island.
    double min_global_ndvi;             // Get the minimum NDVI on all islands.
    Ndvi(Parameters* pars);             // Constructor for the Ndvi class.
    ~Ndvi();                            // Deconstructor.
private:
    void ReadNDVI(Parameters* pars);    // Function to read and store the NDVI data.
                                        // Also sets the simulation length.
};

class Seroprev
{
public:
    int* age_group;             // Age group of the animals, -1 indicates
                                // that the data was not age-stratified.
    int* epi_year_start;        // Epidemiological year the data is collected from.
    int* epi_month;             // Epidemiological month the data is collected from.
    int* n_tested;              // Number of animals tested.
    int* n_positive;            // Number of animals tested which were positive.
    int n_entries;              // The number of entries in the data.
    int* island_id;             // ID number of the island the data entries are from.
    Seroprev();                 // Default constructor of the class.
    ~Seroprev();                // Default destructor of the class.
private:
    void ReadSeroprev();        // Function to read and store the seroprevalence data.
};

class Data
{
public:
    double start_month;         // Starting month of the data.
    double start_year;          // Starting year of the data.
    Ndvi ndvi;                  // NDVI data.
    Seroprev sero;              // Seroprevalence data.
    Data(Parameters* pars);     // Constructor.
private:

};

#endif