#ifndef CONST_H
#define CONST_H 1

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

// Internal unit is 1/h Mpc, km/s, solar mass.
static const double unit_length=   3.085678e24; // Mpc in cm
static const double unit_mass=     1.989e33;    // solar mass in g
static const double unit_velocity= 1.0e5;       // 1 km/s in cm/sec
static const double G_cgs=         6.672e-8;    // Grav const int cgs

static const double unit_time=     unit_length/unit_velocity;
static const double G= G_cgs*unit_mass*(unit_time*unit_time)/
                       (unit_length*unit_length*unit_length);

static const double H0= 100.0;                  // km/s/(1/h Mpc)
static const double rho_crit_0= 3.0*H0*H0/(8.0*M_PI*G);

static const double delta_c= 1.686;
//static const double rho_crit_0= 2.7754e2; // solar_mass/(1/h kpc)^3

#endif
