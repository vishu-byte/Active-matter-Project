#include "../include/Physics.h"
#include "../include/Screen.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_timer.h>
#include <bits/types/time_t.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <time.h>

void state_before_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                             ParSim ::Physics &physics, int steps,
                             int dimension);
void state_after_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                            ParSim ::Physics &physics);

int main() {

  // 1)Creating and initializaing particle system

  /*Parameters*/
  /*Try to stick to S.I units to make sense out of numbers*/
  int Number_of_particles = 2;
  int Number_of_time_steps = 4500;
  int dimension = 500; // meters

  ParSim::ParticleSystem parsym(Number_of_particles);
  ParSim::Physics physics;
  ParSim::Screen screen; // initialize the screen

  ParSim::Particle *const particle =
      parsym.get_particles(); // get access to paticles

  /*Setting physics parameters*/
  physics.parameters[8] = 0.001;           // time step
  physics.parameters[0] = 1000;            // k
  physics.parameters[1] = 2;               // interaction_radius r
  physics.parameters[2] = 0.001;               // mass
  physics.parameters[3] = 2;               // radius
  physics.parameters[4] = 0.8;             // mu
  physics.parameters[5] = 2;               // gamma
  physics.parameters[6] = 0.00000001;      // epsilon1  -- softening length
  physics.parameters[7] = M_PI / 10000000; // epsilon2 -- softening omega

  /*Initial conditions*/
  // particle 1
  particle[0].x = -3;
  particle[0].y = 0;
  particle[0].vx = 6;
  particle[0].vy = 0;
  particle[0].alpha = 0;
  particle[0].omega = -2 * M_PI;
  particle[0].vx_activity = 6;
  particle[0].vy_activity = 0;
  particle[0].omega_activity = 3 * M_PI;

  // particle 2
  particle[1].x = 3;
  particle[1].y = 0;
  particle[1].vx = 0;
  particle[1].vy = 0;
  particle[1].alpha = 0;
  particle[1].omega = 0;
  particle[1].vx_activity = -6;
  particle[1].vy_activity = 0;
  particle[1].omega_activity = 1 * M_PI;

  // 2)Creating a data file for strorage and log-----------

  std::ofstream data_output;
  std::ofstream log;

  data_output.open("data1.xyz");
  log.open("log.txt");

  // Print the state before the simulation in log
  state_before_simulation(log, parsym, physics, Number_of_time_steps,
                          dimension);

  log << "-x-x-x-x-x-Simulation initiated-x-x-x-x-x- " << std::endl;
  std::cout << "-x-x-x-x-x-Simulation initiated-x-x-x-x-x- " << std::endl;

  time_t start = time(&start); // for measuring total runtime

  // 3) Main simulation loop--------------
  for (int step = 0; step < Number_of_time_steps; step++) {

    // Draw current particle system in the window
    screen.draw_particlesystem(parsym);

    // Manipulate particle positions for next iteration.
    physics.evolve_system(parsym, step, log);

    // writing data of this state to file (will be used for rendering the system
    // in ovito)
    data_output << Number_of_particles << std::endl;
    data_output << ' ' << std::endl;

    for (int i = 0; i < parsym.no_of_particles; ++i) {
      data_output << particle[i].x << ' ' << particle[i].y << ' ' << 0 << ' '
                  << cos(particle[i].alpha) << ' ' << sin(particle[i].alpha)
                  << ' ' << 0 << ' ' << particle[i].alpha << ' '
                  << particle[i].vx << ' ' << particle[i].vy << ' '
                  << particle[i].omega << ' ' << std::endl;
    }

    if (step % 100 == 0) {
      std ::cout << "----------Step count: " << step << std::endl;
    }

    log << "----------Step count: " << step << std::endl;
  }

  time_t end = time(&end);

  data_output.close();

  std::cout << "-x-x-x-x-x-Simulation ended-x-x-x-x-x-" << std::endl;
  log << "-x-x-x-x-x-Simulation ended-x-x-x-x-x-" << std::endl;
  log << "Runtime: " << end - start << " seconds" << std::endl;

  /*----------------------------------*/

  // Print the state before the simulation in log
  state_after_simulation(log, parsym, physics);

  return 0;
}

void state_before_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                             ParSim ::Physics &physics, int steps,
                             int dimension) {

  ParSim::Particle *const particle =
      parsym.get_particles(); // get access to paticles
  log << "-------Initial conditions------" << std::endl;

  for (int i = 0; i < parsym.no_of_particles; ++i) {
    log << "Particle: " << i << std::endl;
    log << "x, y = " << particle[i].x << ", " << particle[i].y << std::endl;
    log << "V = " << particle[i].vx << ", " << particle[i].vy << std::endl;
    log << "Omega = " << particle[i].omega << std::endl;
    log << "V0 = " << particle[i].vx_activity << ", " << particle[1].vy_activity
        << std::endl;
    log << "Omega0 = " << particle[i].omega_activity << std::endl;
  }

  log << "-------Parameters and state------" << std::endl;
  log << "Number of particles: " << parsym.no_of_particles << std::endl
      << "Time step: " << physics.parameters[8] << std::endl
      << "Number of time steps: " << steps << std::endl
      << "Dimension: " << dimension << std::endl
      << "k: " << physics.parameters[0] << std::endl
      << "Interaction radius (sigma): " << physics.parameters[1] << std::endl
      << "Radius (r) : " << physics.parameters[3] << std::endl
      << "Mass (m): " << physics.parameters[2] << std::endl
      << "mu: " << physics.parameters[4] << std::endl
      << "gamma: " << physics.parameters[5] << std::endl
      << "epsilon1 " << physics.parameters[6] << std::endl
      << "epsilon2 " << physics.parameters[7] << std::endl;

  log << "Energy-momentum before the collision: " << std::endl;
  log << "Total Energy: "
      << physics.EnergyMomentum(parsym)[0] + physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Translational K.Energy: " << physics.EnergyMomentum(parsym)[0]
      << std::endl;
  log << "Rotational K.Energy: " << physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Momentum: "
      << "(" << physics.EnergyMomentum(parsym)[1] << ", "
      << physics.EnergyMomentum(parsym)[2] << ")" << std::endl;
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
      << "(" << physics.EnergyMomentum(parsym)[1] << ", "
      << physics.EnergyMomentum(parsym)[2] << ")" << std::endl;
}