//
// Created by Andrei Vasilev on 1/13/17.
//

#include "../include/ParticleSystem.h"
#include <SDL2/SDL.h>
#include <cmath>
#include <random>
#include <stdlib.h>

/*Class Particle definitions ------------*/
ParSim::Particle::Particle() { random_initialize(); }

ParSim::Particle::Particle(double x_cor, double y_cor, double v_x, double v_y,
                           double orientation) {
  x = x_cor;
  y = y_cor;
  vx = v_x;
  vy = v_y;
  alpha = orientation;
}

void ParSim::Particle::random_initialize(void) {

  int Number_of_particles = 20;
  int Number_of_time_steps = 1000;
  double phi = 0.60; // area fraction
  double L;
  L = std::sqrt(M_PI * Number_of_particles / phi);
  std::random_device rd;
  std::uniform_real_distribution<double> x_coordinate(-1, 1);
  std::uniform_real_distribution<double> y_coordinate(-1, 1);
  std::uniform_real_distribution<double> vx_dist(-1, 1);
  std::uniform_real_distribution<double> vy_dist(-1, 1);
  std::uniform_real_distribution<double> alpha_dist(0, 1);

  x = L * x_coordinate(rd);
  y = L * y_coordinate(rd); // random distribution

  // Generate random particle speed. Speed is squared causing
  // particle distribution to be exponential instead of linear.
  vx = 3 * pow(vx_dist(rd), 2);
  vy = 3 * pow(vy_dist(rd), 2);

  // Generate random particle orientation (0 to 2pi)
  alpha = alpha_dist(rd);

  // Generate random V0
  vx_activity = 3 * pow(vx_dist(rd), 2);
  vy_activity = 3 * pow(vy_dist(rd), 2);

  // Generatoe random omega
  omega_activity = 2*M_PI*pow(vx_dist(rd), 2);
}

/*Class Particle System definitions----------------*/
ParSim::ParticleSystem::ParticleSystem(int num_of_particles) {
  this->no_of_particles = num_of_particles;
  this->particle_array = new Particle[no_of_particles];
}

ParSim::ParticleSystem::~ParticleSystem() { delete[] particle_array; }

namespace ParSim {
Particle *const ParticleSystem::get_particles() { return particle_array; }
} // namespace ParSim

double ParSim::ParticleSystem::distance(Particle par1, Particle par2) {
  double distance =
      pow(pow((par1.x - par2.x), 2) + pow((par1.y - par2.y), 2), 0.5);
  return distance;
}