#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include "Particle/Particle.h"
#include "global.h"



int WriteData() {
	return DONE;
}



// Function to create a directory

bool createDirectory(const std::string& path) {
	// Create a folder with permissions 0777 (full access for user, group, others)
	int status = mkdir(path.c_str(), 0777);

	if (status == 0) {
		std::cout << "Folder created successfully." << std::endl;
	} else {
		std::cerr << "Error creating folder." << std::endl;
		// You can use perror to print the error message for more details
		perror("mkdir");
	}
	return true;
}



int writeParticle(std::vector<Particle*> &particle, double current_time, int outputNum) {

		const int width = 18;
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

		outputFile << current_time*1e4 << " Myr";
		outputFile << "\n";
    outputFile << std::left << std::setw(width) << "Mass (Msun)"
			<< std::setw(width) << "X (pc)"
			<< std::setw(width) << "Y (pc)"
			<< std::setw(width) << "Z (pc)"
			<< std::setw(width) << "Vx (km/s)"
		 	<< std::setw(width) << "Vy (km/s)" 
			<< std::setw(width) << "Vz (km/s)" << "\n";


    // Write particle data to the file
		for (Particle* ptcl:particle) {
			if (current_time == ptcl->CurrentTimeIrr)
        outputFile  << std::left << std::setw(width) << ptcl->Mass*mass_unit 
                    << std::setw(width) << ptcl->Position[0]*position_unit
                    << std::setw(width) << ptcl->Position[1]*position_unit
                    << std::setw(width) << ptcl->Position[2]*position_unit
                    << std::setw(width) << ptcl->Velocity[0]*velocity_unit/yr*pc/1e5
                    << std::setw(width) << ptcl->Velocity[1]*velocity_unit/yr*pc/1e5 
                    << std::setw(width) << ptcl->Velocity[2]*velocity_unit/yr*pc/1e5 << '\n';
			else
        outputFile  << std::left << std::setw(width) << ptcl->Mass*mass_unit 
                    << std::setw(width) << ptcl->PredPosition[0]*position_unit 
                    << std::setw(width) << ptcl->PredPosition[1]*position_unit 
                    << std::setw(width) << ptcl->PredPosition[2]*position_unit 
                    << std::setw(width) << ptcl->PredVelocity[0]*velocity_unit/yr*pc/1e5 
                    << std::setw(width) << ptcl->PredVelocity[1]*velocity_unit/yr*pc/1e5
                    << std::setw(width) << ptcl->PredVelocity[2]*velocity_unit/yr*pc/1e5 << '\n';
    }

    // Close the file
    outputFile.close();

    std::cout << "Data written to output.txt successfully!" << std::endl;

    return 0;

}
