#ifndef PARTICLE_H
#define PARTICLE_H 1

struct Particle {
  float x[3];
  float z; // debug?
  float vr;
  float ra, dec;
  float w;
};

#endif
