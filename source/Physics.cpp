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

  double k = this->parameters[0];                  // N/m
  double interaction_radius = this->parameters[1]; // m
  double m = this->parameters[2];                  // kg
  double r = this->parameters[3];
  double mu = this->parameters[4];
  double gamma = this->parameters[5];
  double epsilon1 = this->parameters[6];
  double epsilon2 = this->parameters[7];
  double zeta = this->parameters[9]; // zeta
  double fx;
  double fy;

  // Access the particle array. So that notation becomes easier
  Particle *const particle = parsym.get_particles();

  // Loop1: through all particles
  for (int i = 0; i < parsym.no_of_particles; ++i) {

    log << "1st loop -- " << std::endl;

    // Store previous step forces
    particle[i].force_radial_prev[0] = particle[i].force_radial[0];
    particle[i].force_radial_prev[1] = particle[i].force_radial[1];
    particle[i].force_tangential_prev[0] = particle[i].force_tangential[0];
    particle[i].force_tangential_prev[1] = particle[i].force_tangential[1];
    particle[i].torque_prev = particle[i].torque;

    // Reset the current forces
    particle[i].force_radial[0] = 0;
    particle[i].force_radial[1] = 0;
    particle[i].force_tangential[0] = 0;
    particle[i].force_tangential[1] = 0;
    particle[i].torque = 0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added.
    particle[i].force_radial[0] +=
        -2 * zeta * particle[i].vx + particle[i].vx_activity;
    particle[i].force_radial[1] +=
        -2 * zeta * particle[i].vy + particle[i].vy_activity;
    particle[i].torque +=
        -2 * zeta * particle[i].omega + particle[i].omega_activity;

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
        particle[i].force_radial[0] += (interaction_radius - d) *
                                       (particle[i].x - particle[j].x) /
                                       (d + epsilon1);

        particle[i].force_radial[1] += (interaction_radius - d) *
                                       (particle[i].y - particle[j].y) /
                                       (d + epsilon1);

        // tangential friction force
        double N = (interaction_radius -
                    d); // magnitude of radial force used as normal reaction
        double omega_sum = (particle[i].omega + particle[j].omega);

        particle[i].force_tangential[0] +=
            -mu * N * (omega_sum / (abs(omega_sum + epsilon2))) *
            ((particle[i].y - particle[j].y) / (d + epsilon1));

        particle[i].force_tangential[1] +=
            -mu * N * (omega_sum / (abs(omega_sum + epsilon2))) *
            (-(particle[i].x - particle[j].x) / (d + epsilon1));

        // torque on particle

        if (omega_sum != 0) {
          particle[i].torque +=
              -mu * N * (omega_sum / (abs(omega_sum + epsilon2))) * d;
        }

        if (omega_sum == 0) {
          particle[i].torque +=
              -mu * N * (omega_sum / (abs(omega_sum + epsilon2))) * d;
        }

        log << "i: " << i << "j: " << j << "xi: " << parsym.particle_array[i].x
            << "xj: " << parsym.particle_array[j].x
            << "frx_p: " << particle[i].force_radial_prev[0]
            << "fry_p: " << particle[i].force_radial_prev[1]
            << "ftx_p: " << particle[i].force_tangential_prev[0]
            << "fty_p: " << particle[i].force_tangential_prev[1]
            << "tau_p: " << particle[i].torque_prev
            << "frx: " << particle[i].force_radial[0]
            << "fry: " << particle[i].force_radial[1]
            << "ftx: " << particle[i].force_tangential[0]
            << "fty: " << particle[i].force_tangential[1]
            << "tau: " << particle[i].torque << std::endl;
      }
    }
  }
}

void ParSim::Physics::Euler_Integrator(ParSim::Particle &par, int step,
                                       std::ofstream &log) {

  double m = this->parameters[2];
  double time_step = this->parameters[8];
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

  dvx = (Fx)*time_step;
  dvy = (Fy)*time_step;
  log << "dvy: " << dvy << std::endl;
  dw = (Tau)*time_step;

  // update the attributes

  par.x += par.vx * time_step;
  par.y += par.vy * time_step;
  par.vx += dvx;
  par.vy += dvy;

  log << "vy: " << par.vy << std::endl;

  par.alpha += par.omega * time_step;
  par.omega += dw;
}

void ParSim::Physics::Vel_Verlet_Integrator(ParSim::Particle &par, int step,
                                            std::ofstream &log) {

  double m = this->parameters[2];
  double time_step = this->parameters[8];

  // update the attributes

  par.x += (time_step * par.vx) +
           (pow(time_step, 2)) *
               (par.force_radial[0] + par.force_tangential[0]) / (2 * m);
  par.y += (time_step * par.vy) +
           (pow(time_step, 2)) *
               (par.force_radial[1] + par.force_tangential[1]) / (2 * m);

  par.alpha +=
      (time_step * par.omega) + (pow(time_step, 2)) * (par.torque) / (2 * m);

  if (step != 0) {
    par.vx += time_step *
              ((par.force_radial[0] + par.force_tangential[0]) +
               (par.force_radial_prev[0] + par.force_tangential_prev[0])) /
              (2 * m);
    par.vy += time_step *
              ((par.force_radial[1] + par.force_tangential[1]) +
               (par.force_radial_prev[1] + par.force_tangential_prev[1])) /
              (2 * m);
    par.omega += time_step * (par.torque + par.torque_prev) / (2 * m);
  }

  log << "vy: " << par.vy << std::endl;
}

void ParSim::Physics::ERM_Integrator1(ParSim::Particle &par, int step,
                                      std::ofstream &log) {

  double m = this->parameters[2];
  double time_step = this->parameters[8];
  double Fx;
  double Fy;
  double Tau;
  // Total force and torque on this particle
  Fx = par.force_radial[0] + par.force_tangential[0];
  Fy = par.force_radial[1] + par.force_tangential[1]; // -g newton, m = 1
  log << "Ftangential: " << par.force_tangential[1] << std::endl;
  Tau = par.torque;

  // update the attributes upto midpoint

  par.x += par.vx * time_step / 2; // x'
  par.y += par.vy * time_step / 2;
  par.vx += (Fx)*time_step / 2; // v'
  par.vy += (Fy)*time_step / 2;

  par.alpha += par.omega * time_step / 2;
  par.omega += (Tau)*time_step / 2;

  //Error estimation in x and v

}


void ParSim::Physics::ERM_Integrator2(ParSim::Particle &par, int step,
                                      std::ofstream &log) {

  double m = this->parameters[2];
  double time_step = this->parameters[8];
  double Fx;
  double Fy;
  double Tau;
  // Total force and torque on this particle
  Fx = par.force_radial[0] + par.force_tangential[0];
  Fy = par.force_radial[1] + par.force_tangential[1]; // -g newton, m = 1
  log << "Ftangential: " << par.force_tangential[1] << std::endl;
  Tau = par.torque;

  // update the attributes upto midpoint

  par.x += par.vx * time_step; // x'
  par.y += par.vy * time_step;
  par.vx += (Fx)*time_step ; // v'
  par.vy += (Fy)*time_step;

  par.alpha += par.omega * time_step;
  par.omega += (Tau)*time_step ;
}


void ParSim ::Physics::Integrator(ParSim::ParticleSystem &parsym, int step,
                                  std::ofstream &log) {

  for (int i = 0; i < parsym.no_of_particles; ++i) {
    // Vel_Verlet_Integrator(parsym.particle_array[i], step, log);
    ERM_Integrator1(parsym.particle_array[i], step, log);
    //  boundary conditions

    // if (parsym.particle_array[i].x < -1000 ||
    //     parsym.particle_array[i].x > 1000 ||
    //     parsym.particle_array[i].y < -800 || parsym.particle_array[i].y >
    //     800)
    //   parsym.particle_array[i].random_initialize();
  }
}

void ParSim::Physics::evolve_system(ParticleSystem &parsym, int step,
                                    std::ofstream &log) {

  // 1)Force-linking--------
  ParSim::Physics::Force_PP(parsym, log); // links forces on each object
                                          // at this stage
  // 2)Integrating----------
  ParSim::Physics::Integrator(
      parsym, step,
      log); // applies forces on each object as determined above
}

void ParSim::Physics::evolve_system_ERM(ParticleSystem &parsym, int step,
                                        std::ofstream &log) {

  // i) Calculate force (F) from positions and velocities (x,v) --
  Force_PP(parsym, log); // links forces on each object
                                          
  // ii) Update x,v to x', v'  ----
  
  for (int i = 0; i < parsym.no_of_particles; ++i) {
    // Vel_Verlet_Integrator(parsym.particle_array[i], step, log);
    ERM_Integrator1(parsym.particle_array[i], step, log);
    //  boundary conditions

    // if (parsym.particle_array[i].x < -1000 ||
    //     parsym.particle_array[i].x > 1000 ||
    //     parsym.particle_array[i].y < -800 || parsym.particle_array[i].y >
    //     800)
    //   parsym.particle_array[i].random_initialize();
  }

  // iii) Again calculate force (F') from x',v' ------

  Force_PP(parsym, log);

  // error tolerance condition ----

  //  iv) Update x',v' to xnew, vnew --------
  for (int i = 0; i < parsym.no_of_particles; ++i) {
    // Vel_Verlet_Integrator(parsym.particle_array[i], step, log);
    ERM_Integrator2(parsym.particle_array[i], step, log);
    //  boundary conditions

    // if (parsym.particle_array[i].x < -1000 ||
    //     parsym.particle_array[i].x > 1000 ||
    //     parsym.particle_array[i].y < -800 || parsym.particle_array[i].y >
    //     800)
    //   parsym.particle_array[i].random_initialize();
  }

  //v)Change delt t
  //this->parameters[8] = 0.9;
  
}

std::vector<double> ParSim::Physics::EnergyMomentum(ParticleSystem &parsym) {
  std::vector<double> p{0, 0, 0, 0}; // kinetic, rotational, momenta
  double m = this->parameters[2];
  double I = m;
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
