#include "defs.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <sstream>

extern double global_time, EnzoTimeStep;

using TimePoint = std::chrono::high_resolution_clock::time_point;
using Duration  = std::chrono::duration<double>;

struct TimeEntity {
		TimePoint start;
		TimePoint end;
		Duration duration;

		void getDuration() {
			duration += end-start;
		}

		void markStart() {
			start = std::chrono::high_resolution_clock::now();
		}

		void markEnd() {
			end = std::chrono::high_resolution_clock::now();
		}
};

class TimeTracer {
	private:

	public:
		TimeEntity reg;
		TimeEntity reg_sendall;
		TimeEntity reg_gpu;
		TimeEntity reg_cpu1;
		TimeEntity reg_cpu2;
		TimeEntity reg_cpu3;
		TimeEntity reg_cpu4;
		TimeEntity reg_cpu5;

		TimeEntity reg_dt1;
		TimeEntity reg_dt2;
		TimeEntity reg_dt3;

		TimeEntity irr;
		TimeEntity irr_chain;
		TimeEntity irr_force;

	void output() {

		const int width = 25;
		double current_time = global_time*EnzoTimeStep*1e10/1e6;
    std::string fname = "performance.out";


    // Open a file for writing
    std::ofstream outputFile(fname);


    // Check if the file is opened successfully
    if (!outputFile.is_open()) {
			std::cerr << "Error opening the file!" << std::endl;
    }

		outputFile << "This file summarizes the performance of Nbody+ in detail (unit:second).\n"; //
		outputFile << "CurrentTime is " << current_time << " Myr\n"; //

		outputFile << "----------------------------------------------------\n" ;
		outputFile << std::left\
			<< std::setw(width) << "| Regular Routine = " << reg.duration.count() << " |\n";

		outputFile << "----------------------------------------------------\n" ;
		outputFile << std::left\
			<< std::setw(width) << "Regular SendAll" \
			<< std::setw(width) << "Regular GPU" \
			<< std::setw(width) << "Regular CPU1" \
			<< std::setw(width) << "Regular CPU2" \
			<< std::setw(width) << "Regular CPU3" \
			<< std::setw(width) << "Regular CPU4" \
			<< std::setw(width) << "Regular CPU5" \
		 	<< '\n';

		outputFile  << std::left
			<< std::setw(width) << reg_sendall.duration.count() \
			<< std::setw(width) << reg_gpu.duration.count() \
			<< std::setw(width) << reg_cpu1.duration.count() \
			<< std::setw(width) << reg_cpu2.duration.count() \
			<< std::setw(width) << reg_cpu3.duration.count() \
			<< std::setw(width) << reg_cpu4.duration.count() \
			<< std::setw(width) << reg_cpu5.duration.count() \
			<< '\n';

		outputFile << "----------------------------------------------------\n" ;
		outputFile << std::left\
			<< std::setw(width) << "Regular TimeStep1" \
			<< std::setw(width) << "Regular TimeStep2" \
			<< std::setw(width) << "Regular TimeStep3" \
		 	<< '\n';

		outputFile  << std::left
			<< std::setw(width) << reg_dt1.duration.count() \
			<< std::setw(width) << reg_dt2.duration.count() \
			<< std::setw(width) << reg_dt3.duration.count() \
			<< '\n';


		outputFile << "----------------------------------------------------\n" ;
		outputFile << "----------------------------------------------------\n" ;
		outputFile << std::left\
			<< std::setw(width) << "| Irregular Routine = " << irr.duration.count() << " |\n";

		outputFile << "----------------------------------------------------\n" ;
		outputFile << std::left\
			<< std::setw(width) << "Particle Chain" \
			<< std::setw(width) << "Force Calculation" \
		 	<< '\n';

		outputFile  << std::left
			<< std::setw(width) << irr_chain.duration.count() \
			<< std::setw(width) << irr_force.duration.count() \
			<< '\n';

		outputFile << "----------------------------------------------------\n" ;

		outputFile << std::left\
			<< std::setw(width) << "\nComputational Efficiency (Simulation time in Myr/physical time in second)\n";

		outputFile << std::left\
			<< std::setw(width) << "Regular Force" \
			<< std::setw(width) << "Irregular Force" \
		 	<< '\n';

		outputFile  << std::left
			<< std::setw(width) << current_time/reg.duration.count() \
			<< std::setw(width) << current_time/irr.duration.count() \
			<< '\n';
		// Close the file
    outputFile.close();

    //std::cout << "Data written to performance.out successfully!" << std::endl;
	}
};
