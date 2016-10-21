#ifndef PARTICLE_H
#define PARTICLE_H 1

struct Particle {
  float x[3];      // 0,1,2
  float z;         // 3
  float vr;        // 4
  float radec[2];  // 5,6
  float M;         // 7
  float w;         // 8
  float flag;      // 9 
  float rsat, vsat; // 10, 11
};

#endif
