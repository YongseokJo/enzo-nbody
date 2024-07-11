#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include "global.h"

int getLineNumber();
void write_out(std::ofstream& outputFile, const Particle* ptcl);
void write_neighbor(std::ofstream& outputFile, const Particle* ptcl); 
const int NUM_COLUMNS = 7; // Define the number of columns
const int width = 18;



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

    std::string directoryPath = "output";

    // Create the directory or check if it already exists
    if (!createDirectory(directoryPath)) {
        // Handle the error if necessary
        return 1;
    }


    // Now let's save the outputs in a new directory

    // Construct the filename with the timestamp
    std::string fname = directoryPath + "/output_" + std::to_string(outputNum) + ".txt";
    std::string nn_fname = directoryPath + "/neighbor/nn_" + std::to_string(outputNum) + ".txt";

    // Open a file for writing
    std::ofstream output(fname);
    std::ofstream output_nn(nn_fname);


    // Check if the file is opened successfully
    if (!output.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        return 1;
    }

		output << current_time*EnzoTimeStep*1e10/1e6 << " Myr, "; //
		output << global_time*EnzoTimeStep*1e10/1e6 << " Myr"; //
		output << "\n";
		output << outputTime << ", "; //
		output << outputTimeStep << ", "; //
		output << global_time << ""; //
		output << "\n";
    output << std::left 
			<< std::setw(width) << "PID"
			<< std::setw(width) << "Mass (Msun)"
			<< std::setw(width) << "X (pc)"
			<< std::setw(width) << "Y (pc)"
			<< std::setw(width) << "Z (pc)"
			<< std::setw(width) << "Vx (km/s)"
		 	<< std::setw(width) << "Vy (km/s)" 
			<< std::setw(width) << "Vz (km/s)" << "\n";


    // Write particle data to the file
		for (Particle* ptcl:particle) {
			ptcl->predictParticleSecondOrderIrr(current_time);
			if (ptcl->isCMptcl)  {
				ptcl->convertBinaryCoordinatesToCartesian();
				write_out(output, ptcl->BinaryParticleI);
				//write_neighbor(output_nn, ptcl->BinaryParticleI);
				write_out(output, ptcl->BinaryParticleJ);
				//write_neighbor(output_nn, ptcl->BinaryParticleJ);
			}
			else {
				write_out(output, ptcl);
				//write_neighbor(output_nn, ptcl);
			}
    }


    // Close the 
    output.close();
    output_nn.close();

    std::cout << "Data written to output.txt successfully!" << std::endl;

    return 0;

}


void write_out(std::ofstream& outputFile, const Particle* ptcl) {
        outputFile  << std::left
										<< std::setw(width) << ptcl->PID
										<< std::setw(width) << ptcl->Mass*mass_unit
                    << std::setw(width) << ptcl->PredPosition[0]*position_unit
                    << std::setw(width) << ptcl->PredPosition[1]*position_unit
                    << std::setw(width) << ptcl->PredPosition[2]*position_unit
                    << std::setw(width) << ptcl->PredVelocity[0]*velocity_unit/yr*pc/1e5
                    << std::setw(width) << ptcl->PredVelocity[1]*velocity_unit/yr*pc/1e5
                    << std::setw(width) << ptcl->PredVelocity[2]*velocity_unit/yr*pc/1e5
										<< std::endl;


}


void write_neighbor(std::ofstream& outputFile, const Particle* ptcl) {
	outputFile  << std::left\
			<< std::setw(width) << ptcl->PID << " = " ;
	for (Particle* nn:ptcl->ACList) {
			outputFile << nn->PID << ", ";
	
	}
	outputFile << "\n";

}

#ifdef time_trace
void output_time_trace() {

}
#endif
