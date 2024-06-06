#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sys/stat.h>
#include "Particle.h"
#include "global.h"

// Function to create a directory

bool createDirectory(const std::string& path) {
    #ifdef _WIN32
        if (_mkdir(path.c_str()) == 0)
    #else
        if (mkdir(path.c_str(), 0777) == 0)
    #endif
        {
            std::cout << "Directory created successfully!" << std::endl;
            return true;
        } else {
            // Check if the directory already exists
            struct stat info;
            if (stat(path.c_str(), &info) == 0 && (info.st_mode & S_IFDIR) != 0) {
                std::cout << "Directory already exists." << std::endl;
                return true;
            }

            std::cerr << "Error creating the directory!" << std::endl;
            return false;
        }
}



int output(std::vector<Particle*> &particle, int outputNum) {

    // First let's create a new directory
    // only create a new directory if the directory doesn't exist

    std::string directoryPath = "output";

    // Create the directory or check if it already exists
    if (!createDirectory(directoryPath)) {
        // Handle the error if necessary
        return 1;
    }


    // Now let's save the outputs in a new directory

    // Construct the filename with the timestamp
    std::string filename = directoryPath + "/snapshot_" +std::to_string(outputNum) + ".txt";

    // Open a file for writing
    std::ofstream outputFile(filename);

    // Check if the file is opened successfully
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        return 1;
    }

    outputFile << "Mass\tX\tY\tZ\tVx\tVy\tVz\n";


    // Write particle data to the file
    for (int i = 0; i < NNB; ++i) {
        outputFile  << particle[i]->Mass << '\t'
                    << particle[i]->Position[0] << '\t'
                    << particle[i]->Position[1] << '\t'
                    << particle[i]->Position[2] << '\t'
                    << particle[i]->Velocity[0] << '\t'
                    << particle[i]->Velocity[1] << '\t'
                    << particle[i]->Velocity[2] << '\n';
    }

    // Close the file
    outputFile.close();

    std::cout << "Data written to output.txt successfully!" << std::endl;

    return 0;

}