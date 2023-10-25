#include "../include/Physics.h"
#include <bits/types/time_t.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <time.h>

void state_before_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                             ParSim ::Physics &physics, int steps,
                             double dimension, double phi);
void state_after_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                            ParSim ::Physics &physics);

int main() {

  // 1)Creating and initializaing particle system

  /*Parameters*/
  /*Try to stick to S.I units to make sense out of numbers*/
  int Number_of_particles = 2;
  int Number_of_time_steps = 4000;
  double phi = 0.60; // area fraction
  double L;
  L = std::sqrt(M_PI * Number_of_particles / phi);
  L = 18;
  ParSim::ParticleSystem parsym(Number_of_particles, phi, L);
  ParSim::Physics physics;

  ParSim::Particle *const particle =
      parsym.get_particles(); // get access to paticles

  /*Setting physics parameters*/
  physics.parameters[8] = 0.001;           // time step
  physics.parameters[0] = 1500;            // k
  physics.parameters[1] = 1;               // interaction_radius sigma
  physics.parameters[2] = 1;               // mass
  physics.parameters[3] = 1;               // radius
  physics.parameters[4] = 0.8;             // mu
  physics.parameters[5] = 0.0;             // gamma
  physics.parameters[6] = 0.00000001;      // epsilon1  -- softening length
  physics.parameters[7] = M_PI / 10000000; // epsilon2 -- softening omega
  physics.parameters[9] = 0.5 * physics.parameters[5] /
                          pow((physics.parameters[2] * physics.parameters[0]),
                              0.5); // zeta

  physics.parameters[10] = 50; // eta      --increase judiciuosly, it should not overpower k
  physics.parameters[11] = 0;   // Diffusion constant

  /*Initial conditions*/
  // particle 1
  particle[0].x = 0;
  particle[0].y = 5;
  particle[0].vx = 6;
  particle[0].vy = 0;
  particle[0].alpha = 0;
  particle[0].omega = 0;
  particle[0].vx_activity = 0;
  particle[0].vy_activity = 0;
  particle[0].omega_activity = 0 * M_PI;

  // particle 2
  particle[1].x = 0;
  particle[1].y = 0;
  particle[1].vx = 0;
  particle[1].vy = 0;
  particle[1].alpha = 0;
  particle[1].omega = 0;
  particle[1].vx_activity = 0;
  particle[1].vy_activity = 0;
  particle[1].omega_activity = 0 * M_PI;

  // 2)Creating a data file for storage and log-----------

  std::ofstream data_output;
  std::ofstream log;

  data_output.open("data1.xyz");
  log.open("log.txt");

  // Print the state before the simulation in log
  state_before_simulation(log, parsym, physics, Number_of_time_steps, L, phi);

  log << "-x-x-x-x-x-Simulation initiated-x-x-x-x-x- " << std::endl;
  std::cout << "-x-x-x-x-x-Simulation initiated-x-x-x-x-x- " << std::endl;

  time_t start = time(&start); // for measuring total runtime

  // 3) Main simulation loop--------------
  for (int step = 0; step < Number_of_time_steps; step++) {

    // writing data of this state to file (will be used for rendering the system
    // in ovito)

    data_output << Number_of_particles << std::endl;
    data_output << "Lattice="
                << "\"10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 0.0\"" << std::endl;
    // first store current configuration
    for (int i = 0; i < parsym.no_of_particles; ++i) {
      data_output << particle[i].x << ' ' << particle[i].y << ' ' << 0 << ' '
                  << cos(particle[i].alpha) << ' ' << sin(particle[i].alpha)
                  << ' ' << 0 << ' ' << particle[i].alpha << ' '
                  << particle[i].vx << ' ' << particle[i].vy << ' '
                  << particle[i].omega << ' ' << std::endl;
    }

    if (step % 50 == 0) {
      std ::cout << "----------Step count: " << step << std::endl;
    }

    log << "----------Step count: " << step << std::endl;

    // Manipulate particle positions for next iteration.
    physics.evolve_system_ERM(parsym, step, log);
  }

  time_t end = time(&end);

  data_output.close();

  std::cout << "-x-x-x-x-x-Simulation ended-x-x-x-x-x-" << std::endl;
  log << "-x-x-x-x-x-Simulation ended-x-x-x-x-x-" << std::endl;
  log << "Runtime: " << end - start << " seconds" << std::endl;

  /*----------------------------------*/

  // Print the state before the simulation in log
  state_after_simulation(log, parsym, physics);

  // log.close();
  // logv.close();
  // logx.close();

  return 0;
}

void state_before_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                             ParSim ::Physics &physics, int steps,
                             double dimension, double phi) {

  ParSim::Particle *const particle =
      parsym.get_particles(); // get access to paticles

  log << "-------Parameters and state------" << std::endl;
  log << "Number of particles: " << parsym.no_of_particles << std::endl
      << "Time step: " << physics.parameters[8] << std::endl
      << "Number of time steps: " << steps << std::endl
      << "Phi: " << phi << std::endl
      << "Dimension: " << dimension << std::endl
      << "k: " << physics.parameters[0] << std::endl
      << "Interaction radius (sigma): " << physics.parameters[1] << std::endl
      << "Radius (r) : " << physics.parameters[3] << std::endl
      << "Mass (m): " << physics.parameters[2] << std::endl
      << "mu: " << physics.parameters[4] << std::endl
      << "gamma: " << physics.parameters[5] << std::endl
      << "zeta: " << physics.parameters[9] << std::endl
      << "eta: " << physics.parameters[10] << std::endl
      << "D: " << physics.parameters[11] << std::endl
      << "epsilon1 " << physics.parameters[6] << std::endl
      << "epsilon2 " << physics.parameters[7] << std::endl
      << "seed " << physics.parameters[7] << std::endl;

  log << "Energy-momentum before the collision: " << std::endl;
  log << "Total Energy: "
      << physics.EnergyMomentum(parsym)[0] + physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Translational K.Energy: " << physics.EnergyMomentum(parsym)[0]
      << std::endl;
  log << "Rotational K.Energy: " << physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Momentum: "
      << "(" << physics.EnergyMomentum(parsym)[2] << ", "
      << physics.EnergyMomentum(parsym)[3] << ")" << std::endl;

  log << "-------Initial conditions------" << std::endl;

  for (int i = 0; i < parsym.no_of_particles; ++i) {
    log << "Particle: " << i << " ---------" << std::endl;
    log << "x, y = " << particle[i].x << ", " << particle[i].y << std::endl;
    log << "V = " << particle[i].vx << ", " << particle[i].vy << std::endl;
    log << "Omega = " << particle[i].omega << std::endl;
    log << "V0 = " << particle[i].vx_activity << ", " << particle[1].vy_activity
        << std::endl;
    log << "Omega0 = " << particle[i].omega_activity << std::endl;
  }
};

void state_after_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                            ParSim ::Physics &physics) {
  log << "Energy-momentum After the collision: " << std::endl;
  log << "Total Energy: "
      << physics.EnergyMomentum(parsym)[0] + physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Translational K.Energy: " << physics.EnergyMomentum(parsym)[0]
      << std::endl;
  log << "Rotational K.Energy: " << physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Momentum: "
      << "(" << physics.EnergyMomentum(parsym)[2] << ", "
      << physics.EnergyMomentum(parsym)[3] << ")" << std::endl;
}
