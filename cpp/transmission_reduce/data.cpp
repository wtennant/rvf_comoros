// data.cpp: Contains the constructor and destructor of the data classes.
#include <fstream>      // For reading in the data from file.
#include <iostream>     // Input-output to console.
#include "config.h"     // Configuration for the RVF model.
#include "data.h"       // Definition of the data classes.

// Function to read in the NDVI data.
void Ndvi::ReadNDVI(Parameters* pars)
{
    // Create a file stream.
    std::fstream file;
    
    // Open the file to read in.
    file.open(static_cast<std::string>(RVF_IDIR) + "ndvi.csv",
              std::fstream::in);

    // Check that the file opened.
    if (file.is_open())
    {
        // Define the contents of a line.
        std::string line;

        // Throw away the first line of the data (the header).
        getline(file, line);

        // Get the number of lines in the file.
        int n_lines = 0;
        while (getline(file, line))
        {
            ++n_lines;
        }

        // Set the length of the simulation.
        // There are n_island times as many lines as simulation steps.
        pars->n_steps = n_lines / n_island;

        // For each island, allocate space to store the data based on the number
        // of lines in the file. The data is in months, so multiply by four to
        // get the data in weeks.
        for (int i = 0; i < n_island; ++i)
        {
            ndvi[i] = new double[pars->n_steps];
        }

        // Reset the position of the file stream to the first line.
        file.clear();
        file.seekg(0, std::ios::beg);

        // Throw away the first line of the data (the header) again.
        getline(file, line);

        // Keep track of the number of entries scanned, and the
        // number of entries read in to the data.
        int entries = 0;
        int* read_entries = new int[n_island]{0};

        // Define the number of columns in the data.
        int n_cols = 6;

        // Define the separator of entries.
        char sep = ',';

        // Variable to keep track of the island ID number.
        int island_id = 0;

        // For each line, extract the final entry in the data.
        // The final entry is the NDVI data.
        while(getline(file, line, sep))
        {       
            // Increase the counter for the number of entries that have been read.
            ++entries;

            // If this value is divisible by four, then the current line contains
            // the ndvi data.
            // On the first pass, record the month and year that the NDVI data
            // starts at.
            if (entries % n_cols == 0)
            {
                // Write the entry to the NDVI time-series for each
                // epidemiological week.
                ndvi[island_id][read_entries[island_id]] = std::stod(line);

                // Store the minimum NDVI value if necessary.
                if (ndvi[island_id][read_entries[island_id]] < min_local_ndvi[island_id])
                {
                    min_local_ndvi[island_id] = ndvi[island_id][read_entries[island_id]];

                    // Also check if this is a new global minimum.
                    if (min_local_ndvi[island_id] < min_global_ndvi)
                    {
                        min_global_ndvi = min_local_ndvi[island_id];
                    }
                }

                // Increase the number of entries that have been read.
                ++read_entries[island_id];
            }
            else if (entries == 4)
            {
                // Get the starting month of the data.
                start_month = std::stoi(line);
            }
            else if (entries == 5)
            {
                // Get the starting year of the data.
                start_year = std::stoi(line);
            }
            else if ((entries - 1) % n_cols == 0)
            {
                // Get the island ID number for the row.
                island_id = std::stoi(line);
            }

            // Define the next separator character.
            // If about to read the final entry of a line, change the separator
            // to a new line.
            sep = (entries + 1) % n_cols == 0 ? '\n' : ',';
        }

        // Close the file.
        file.close();
    }
    else
    {
        // Write an error messsage and terminate the execution of the program.
        std::cout << "Error: NDVI file could not be read." << std::endl;
        exit(1);
    }
}

// Function to read in the seroprevalence data.
void Seroprev::ReadSeroprev()
{
    // Create a file stream.
    std::fstream file;

    // Open the file to read in.
    file.open(static_cast<std::string>(RVF_IDIR) + "seroprev.csv",
              std::fstream::in);

    // Check that the file opened.
    if (file.is_open())
    {
        // Define the contents of a line.
        std::string line;

        // Throw away the first line of the data (the header).
        getline(file, line);

        // Get the number of lines in the file.
        int n_lines = 0;
        while (getline(file, line))
        {
            ++n_lines;
        }

        // Allocate space to store the data based on the number of lines in the file.
        age_group = new int[n_lines];
        epi_year_start = new int[n_lines];
        n_tested = new int[n_lines];
        n_positive = new int[n_lines];
        island_id = new int[n_lines];
        epi_month = new int[n_lines];
        n_entries = n_lines;

        // Reset the position of the file stream to the first line.
        file.clear();
        file.seekg(0, std::ios::beg);

        // Throw away the first line of the data (the header) again.
        getline(file, line);

        // Define the entry number and line number.
        int entry = 0, i_line = 0;

        // Define the number of columns in the data.
        int n_cols = 7;

        // Define the separator of entries.
        char sep = ',';

        // For each line, extract the final entry in the data.
        // The final entry is the NDVI data.
        while (getline(file, line, sep))
        {
            // Store the data depending on the column that is currently being read.
            if (entry % n_cols == 0)        // First column of the data.
            {
                island_id[i_line] = std::stoi(line);
            }
            else if (entry % n_cols == 1)   // Second column of the data.
            {
                n_positive[i_line] = std::stoi(line);
            }
            else if (entry % n_cols == 2)   // Third column of the data.
            {
                n_tested[i_line] = std::stoi(line);
            }
            else if (entry % n_cols == 3)   // Fourth column of the data.
            {
                age_group[i_line] = std::stoi(line);
            }
            else if (entry % n_cols == 4)   // Fifth column of the data.
            {
                epi_year_start[i_line] = std::stoi(line);
            }
            else if (entry % n_cols == 6)   // Seventh column.
            {
                epi_month[i_line] = std::stoi(line);
            }

            // Define the next separator character.
            // If about to read the final entry of a line, change the separator
            // to a new line.
            sep = (entry + 2) % n_cols == 0 ? '\n' : ',';

            // Increase the line counter if on the final entry of the line.
            i_line = (entry + 1) % n_cols == 0 ? i_line + 1 : i_line;

            // Move on to the next entry.
            ++entry;
        }

        // Close the file.
        file.close();
    }
    else
    {
        // Write an error messsage and terminate the execution of the program.
        std::cout << "Error: seroprevalence data file could not be read." << std::endl;
        exit(3);
    }
}

// Constructor for the seroprevalence class.
Ndvi::Ndvi(Parameters* pars)
{
    // Store the number of islands.
    n_island = pars->n_island;

    // Set up space to store the NDVI data for each island.
    ndvi = new double*[n_island];

    // Set up space to store the minimum NDVI across all islands and for each island.
    min_global_ndvi = 1.5;
    min_local_ndvi = new double[n_island];

    // ... and initialise.
    for (int i = 0; i < n_island; ++i)
    {
        min_local_ndvi[i] = 1.5;
    }

    // Set up space to store the NDVI data and read in the NDVI data.
    // Also set the length of the simulation.
    ReadNDVI(pars);
}

// Deconstructor for the seroprevalence class.
Ndvi::~Ndvi()
{
    // De-allocate space for the seroprevalence data.
    delete[] min_local_ndvi;
    for (int i = 0; i < n_island; ++i)
    {
        delete[] ndvi[i];
    }
    delete[] ndvi;
}

// Constructor for the seroprevalence class.
Seroprev::Seroprev()
{
    // Allocate, read and store the seroprevalence data.
    ReadSeroprev();
}

// Deconstructor for the seroprevalence class.
Seroprev::~Seroprev()
{
    // De-allocate space for the seroprevalence data.
    delete[] island_id;
    delete[] n_positive;
    delete[] n_tested;
    delete[] epi_year_start;
    delete[] age_group;
}

// Constructor for the data class.
Data::Data(Parameters* pars) : ndvi(pars), sero()
{
    // Set the starting year and month equal to the first
    // entry in the NDVI data.
    start_month = ndvi.start_month;
    start_year = ndvi.start_year;
}