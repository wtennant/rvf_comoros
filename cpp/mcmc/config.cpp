// config.cpp: Contains functions related to the configuration of the model
// and fitting procedure.
#include <fstream>      // Writing to file.
#include <iostream>     // Input-output stream.
#include "config.h"     // Configuration definitions.

// Function to output the configuration to file.
void WriteConfig()
{
// Delcare the file stream.
    std::fstream file;

    // Open the file for writing.
    file.open(static_cast<std::string>(RVF_ODIR) + "config.csv", std::fstream::out);

    // Check that the file opened.
    if (file.is_open())
    {
        // If the file could be opened, write the header.
        file << "CONFIG_NAME,CONFIG_VALUE";
        
        // For each configuration, record its setting.
        file << "\nRVF_CHAINS," << RVF_CHAINS;
        file << "\nRVF_MIN_NDVI," << RVF_MIN_NDVI;
        file << "\nRVF_TRANSMISSION," << RVF_TRANSMISSION;
        file << "\nRVF_DIFF_SCALE," << RVF_DIFF_SCALE;
        file << "\nRVF_DIFF_NDVI," << RVF_DIFF_NDVI;

        // Close the file once completed.
        file.close();
    }
    else
    {
        std::cout << "Configuration file could not be opened for writing." << std::endl;
        exit(1);
    }    
}