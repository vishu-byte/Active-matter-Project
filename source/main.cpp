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

int main() {

  // 1)Creating and initializaing particle system

  /*Parameters*/
  /*Try to stick to S.I units to make sense out of numbers*/
  int Number_of_particles = 2;
  double Time_step = 0.001; // seconds
  int Number_of_time_steps = 3500;
  int dimension = 500; // meters

  ParSim::ParticleSystem parsym(Number_of_particles);
  ParSim::Physics physics;
  ParSim::Screen screen; // initialize the screen

  ParSim::Particle *const particle =
      parsym.get_particles(); // get access to paticles

  /*Setting physics parameters*/
  physics.force_params[0] = 1000;        // k
  physics.force_params[1] = 2;           // interaction_radius r
  physics.force_params[2] = 1;           // mass
  physics.force_params[3] = 2;           // radius
  physics.force_params[4] = 0.8;         // mu
  physics.force_params[5] = 0.5;         // gamma
  physics.force_params[6] = M_PI / 2;    // omega_intrinsic
  physics.force_params[7] = M_PI / 1000; // omega_tolerance

  /*Initial conditions*/
  // particle 1
  particle[0].x = -5;
  particle[0].y = 0;
  particle[0].vx = 3;
  particle[0].vy = 0;
  particle[0].alpha = 0;
  particle[0].omega = +2 * M_PI;
  particle[0].vx_activity = 1;
  particle[0].vy_activity = 2;
  particle[0].omega_activity = M_PI;
  


  // particle 2
  particle[1].x = 5;
  particle[1].y = 0;
  particle[1].vx = -3;
  particle[1].vy = 0;
  particle[1].alpha = 0;
  particle[1].omega = 0 * M_PI;
  particle[1].vx_activity = 1;
  particle[1].vy_activity = 2;
  particle[1].omega_activity = M_PI;
  

  // 2)Creating a data file for strorage and log-----------

  std::ofstream data_output;
  std::ofstream log;

  data_output.open("data1.xyz");
  log.open("log.txt");

  log << "-------Parameters and state------" << std::endl;
  log << "Number of particles: " << Number_of_particles << std::endl
      << "Time step: " << Time_step << std::endl
      << "Number of time steps: " << Number_of_time_steps << std::endl
      << "Dimension: " << dimension << std::endl
      << "k: " << physics.force_params[0] << std::endl
      << "Interaction radius (sigma): " << physics.force_params[1] << std::endl
      << "Radius (r) : " << physics.force_params[2] << std::endl
      << "Mass (m): " << physics.force_params[3] << std::endl
      << "mu: " << physics.force_params[4] << std::endl
      << "gamma: " << physics.force_params[5] << std::endl
      << "Intrinsic omega: " << physics.force_params[6] << std::endl
      << "Omega_tolerance: " << physics.force_params[7] << std::endl;

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
  log << "-x-x-x-x-x-Simulation initiated-x-x-x-x-x- " << std::endl;

  std::cout << "-x-x-x-x-x-Simulation initiated-x-x-x-x-x- " << std::endl;

  time_t start = time(&start); // for measuring total runtime

  // 3) Main simulation loop--------------
  for (int step = 0; step < Number_of_time_steps; step++) {

    // Draw current particle system in the window
    screen.draw_particlesystem(parsym);
    // SDL_Delay(10000);

    // Manipulate particle positions for next iteration.
    physics.evolve_system(parsym, Time_step, log);

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

  return 0;
}
