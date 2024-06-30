#include <mpi.h>
#include <iostream>
#include "global.h"
#include "defs.h"


void CalculateSingleAcceleration(Particle *ptcl1, Particle *ptcl2, double (&a)[3], double (&adot)[3]);

int RemoveParticles(std::vector<Particle*> &particle, std::vector<Particle*> &EscapeList) {

	/*
	// Force correction might not needed but neighborhood correction may needed?
	for (Particle* ptcl1:particle) {
		int i=0;
		for (Particle* ptcl2:EscapeList) {
			CalculateSingleAcceleration(ptcl1, ptcl2)
		}
	}
	*/

	// Neighborhood correction
	for (Particle* elem2:EscapeList) {
		int i=0;
		for (Particle* elem1:particle) {
			if (elem1 == elem2) {
				particle.erase(particle.begin() + i);
				break;
			}
		}
	}



}
