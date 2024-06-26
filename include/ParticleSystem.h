//
// Created by Vishu Saini on 31/08/23
//

#ifndef PARTICLE_SIMULATION_PARTICLE_SYSTEM_H
#define PARTICLE_SIMULATION_PARTICLE_SYSTEM_H

#include <cmath>
#include <math.h>
#include <random>
#include <vector>

namespace ParSim { // for particle simulation

class Particle {

public:
  // cartesian coordinates
  double x;
  double y;
  double vx;
  double vy;
  double alpha;
  double omega;
  double radius;
  double vx_activity;
  double vy_activity;
  double omega_activity;
  std::vector<double> force_radial{0, 0};     // radial force: fx and fy
  std::vector<double> force_tangential{0, 0}; // radial force: fx and fy
  double torque{0};                           // torque

  // For storing forces of one step before
  std::vector<double> position_prev{0, 0}; // (x, y) 
  std::vector<double> velocity_prev{0, 0}; // (vx, vy)
  double alpha_prev;
  double omega_prev;

  std::vector<double> force_radial_prev{0, 0};     // radial force: fx and fy
  std::vector<double> force_tangential_prev{0, 0}; // radial force: fx and fy
  double torque_prev{0};                           // torque
  // std:: vector <double> IC{0,0,0,0,0};

  Particle(); // default constructor
  Particle(double x_cor, double y_cor, double speed, double direction,
           double orientation); // parameterized constructor
  virtual ~Particle(){};        // virtual destructor
  void random_initialize(void); // randomly initializes particle
};

class ParticleSystem {
public:
  int no_of_particles;
  Particle *particle_array{nullptr}; // creating particle array on heap

  ParticleSystem(int);             // parameterized constructor
  virtual ~ParticleSystem();       // destructor
  Particle *const get_particles(); // constant pointer, can not change address
                                   // of memory block to which it points
  double distance(Particle par1, Particle par2);
};

} // namespace ParSim

#endif // PARTICLE_SIMULATION_PARTICLE_AND_SYSTEM_H
