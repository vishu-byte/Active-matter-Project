//
// Created by Andrei Vasilev on 1/13/17.
//

#include "../include/Physics.h"
#include "../include/ParticleSystem.h"
#include <SDL2/SDL.h>
#include <cmath>
#include <fstream>
#include <math.h>
#include <random>
#include <stdlib.h>
#include <vector>

// void ParSim::Physics::default_physics(){
//  move(elapsed);
// }

ParSim::Physics::Physics() { int sample = 0; }

/*Force linker + integrators-- */

void ParSim::Physics::Force_PP(ParSim::ParticleSystem &parsym,
                               std::ofstream &log) {
  // parameters of force

  double k = this->force_params[0];                  // N/m
  double interaction_radius = this->force_params[1]; // m
  double m = this->force_params[2];                  // kg
  double r = this->force_params[3];
  double mu = this->force_params[4];
  double gamma = this->force_params[5];
  double epsilon1 = this->force_params[6];
  double epsilon2 = this->force_params[6];

  double fx;
  double fy;

  // Access the particle array. So that notation becomes easier
  Particle *const particle = parsym.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < parsym.no_of_particles; ++i) {

    log << "1st loop -- " << std::endl;
    particle[i].force_radial[0] = 0; // reset the radial force calculator
    particle[i].force_radial[1] = 0;
    particle[i].force_tangential[0] =
        0; // reset the tangential force calculator
    particle[i].force_tangential[1] = 0;
    particle[i].torque = 0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added.
    particle[i].force_radial[0] +=
        -gamma * particle[i].vx + particle[i].vx_activity;
    particle[i].force_radial[1] +=
        -gamma * particle[i].vy + particle[i].vy_activity;
    particle[i].torque +=
        -gamma * particle[i].omega + particle[i].omega_activity;

    // Knary force calculation --- Loop2: through all particles
    for (int j = 0; j < parsym.no_of_particles; ++j) {
      if (j == i) { // no self coupling
        continue;
      }
      log << "2nd loop -- " << std::endl;

      double d = parsym.distance(particle[i], particle[j]);

      // U

      if (d <= interaction_radius) {
        // radial interaction force
        particle[i].force_radial[0] += k * (interaction_radius - d) *
                                       (particle[i].x - particle[j].x) /
                                       (d + epsilon1);

        particle[i].force_radial[1] += k * (interaction_radius - d) *
                                       (particle[i].y - particle[j].y) /
                                       (d + epsilon1);

        // tangential friction force
        double N = k * (interaction_radius -
                        d); // magnitude of radial force used as normal reaction
        double omega_sum = (particle[i].omega + particle[j].omega);

        particle[i].force_tangential[0] +=
            -mu * N * (omega_sum / (abs(omega_sum + epsilon2))) *
            ((particle[i].y - particle[j].y) / (d + epsilon1));

        particle[i].force_tangential[1] +=
            -mu * N * (omega_sum / (abs(omega_sum + epsilon2))) *
            (-(particle[i].x - particle[j].x) / (d + epsilon1));

        // torque on particle

        particle[i].torque +=
            -mu * N * (omega_sum / (abs(omega_sum + epsilon2))) * d;

        log << "i: " << i << "j: " << j << "xi: " << parsym.particle_array[i].x
            << "xj: " << parsym.particle_array[j].x
            << "frx: " << particle[i].force_radial[0]
            << "fry: " << particle[i].force_radial[1]
            << "ftx: " << particle[i].force_tangential[0]
            << "fty: " << particle[i].force_tangential[1]
            << "tau: " << particle[i].torque << std::endl;
      }
    }
  }
}

void ParSim::Physics::Euler_Integrator(ParSim::Particle &par, double time_step,
                                       std::ofstream &log) {

  double m = this->force_params[2];
  double Fx;
  double Fy;
  double Tau;
  double dvx;
  double dvy;
  double dw;
  // Total force and torque on this particle
  Fx = par.force_radial[0] + par.force_tangential[0];
  Fy = par.force_radial[1] + par.force_tangential[1]; // -g newton, m = 1
  log << "Ftangential: " << par.force_tangential[1] << std::endl;
  Tau = par.torque;

  dvx = (Fx / m) * time_step;
  dvy = (Fy / m) * time_step;
  log << "dvy: " << dvy << std::endl;
  dw = (Tau / m) * time_step;

  // update the attributes

  par.x += par.vx * time_step;
  par.y += par.vy * time_step;
  par.vx += dvx;
  par.vy += dvy;

  log << "vy: " << par.vy << std::endl;

  par.alpha += par.omega * time_step;
  par.omega += dw;
}

void ParSim ::Physics::Integrator(ParSim::ParticleSystem &parsym,
                                  double time_step, std::ofstream &log) {

  // parameters
  double m = this->force_params[2];
  for (int i = 0; i < parsym.no_of_particles; ++i) {
    Euler_Integrator(parsym.particle_array[i], time_step, log);

    // boundary conditions

    // if (parsym.particle_array[i].x < -1000 ||
    //     parsym.particle_array[i].x > 1000 ||
    //     parsym.particle_array[i].y < -800 || parsym.particle_array[i].y >
    //     800)
    //   parsym.particle_array[i].random_initialize();
  }
}

void ParSim::Physics::evolve_system(ParticleSystem &parsym, double time_step,
                                    std::ofstream &log) {

  // 1)Force-linking--------
  ParSim::Physics::Force_PP(parsym, log); // links forces on each object
                                          // at this stage
  // 2)Integrating----------
  ParSim::Physics::Integrator(
      parsym, time_step,
      log); // applies forces on each object as determined above
}

std::vector<double> ParSim::Physics::EnergyMomentum(ParticleSystem &parsym) {
  std::vector<double> p{0, 0, 0, 0}; // kinetic, rotational, momenta
  double m = this->force_params[2];
  double I = 1;
  for (int i = 0; i < parsym.no_of_particles; ++i) {
    p[0] += 0.5 * m *
            (pow(parsym.particle_array[i].vx, 2) +
             pow(parsym.particle_array[i].vy, 2)); // KE
    p[1] +=
        0.5 * I * (pow(parsym.particle_array[i].omega, 2)); // Rotational K.E
    p[1] += m * parsym.particle_array[i].vx;                // px
    p[2] += m * parsym.particle_array[i].vy;                // py
  }
  return p;
}
