#define no_well
#ifdef well
#include <iostream>
#include <vector>
#include <algorithm>
#include "../global.h"

class ParticleSystem {
	public:
		void addParticle(const Particle& particle) {
			particles.push_back(particle);
		}

		void removeInactiveParticles() {
			particles.erase(
					std::remove_if(particles.begin(), particles.end(),
						[](const Particle& p) { return !p.active; }),
					particles.end());
		}

		void updateParticles(float deltaTime) {
			for (auto& particle : particles) {
				if (particle.active) {
					particle.x += particle.vx * deltaTime;
					particle.y += particle.vy * deltaTime;
					particle.z += particle.vz * deltaTime;

					if (particle.x < 0 || particle.y < 0 || particle.z < 0) {
						particle.active = false;
					}
				}
			}
			removeInactiveParticles();
		}

		void printParticles() const {
			for (const auto& particle : particles) {
				std::cout << "Particle at (" << particle.x << ", " << particle.y << ", " << particle.z << ")\n";
			}
		}

	private:
		std::vector<Particle> particles;
};


#endif
