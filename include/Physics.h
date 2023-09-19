//
// Created by Vishu Saini on 31/08/23
//

#ifndef PARTICLE_SIMULATION_PHYSICS_H
#define PARTICLE_SIMULATION_PHYSICS_H

#include "ParticleSystem.h"
#include <math.h>
#include <random>

#include <fstream>
#include <iostream> //for writing data
#include <vector>

namespace ParSim { // for particle simulation

class Physics {
public:
  // class responsible for handling all physics behind the simulation
  std::vector<double> force_params{0, 0, 0, 0, 0, 0, 0, 0};

  Physics(); // constructor

  virtual ~Physics(){}; // destructor
public:
  /*Force linker + integrators-- */
  void Force_PP(ParticleSystem &, std::ofstream &); // Main force linker
  void Integrator(ParticleSystem &, double, std::ofstream &); // Main integrator
  void Euler_Integrator(Particle &, double, std::ofstream &);
  void RK4_Integrator(Particle &, int);

  /*Conserved quantities*/
  std::vector<double> EnergyMomentum(ParticleSystem &);

  void evolve_system(ParticleSystem &, double,
                     std::ofstream &file); // takes a particle system and moves
                                           // it forward in time
  void write_to_file(std::ofstream &file, int index, Particle &, int time_step);
};

} // namespace ParSim

#endif // PARTICLE_SIMULATION_PHYSICS_H
