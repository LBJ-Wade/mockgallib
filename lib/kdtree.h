#ifndef KDTREE_H
#define KDTREE_H 1

#include <cstdlib>
#include "particle.h"

typedef int index_t;

struct KDTree {
  float left[3], right[3];
  KDTree* subtree[2];
  Particle* particles[2];
};


size_t kdtree_construct(KDTree* const tree, Particle* const particle, const index_t np, const int quota);

#endif
