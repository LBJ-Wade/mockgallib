#ifndef LIGHTCONE_H
#define LIGHTCONE_H 1

#include <vector>
#include <deque>
#include "halo.h"

class LightCone : public std::vector<Halo> {

};

class LightCones : public std::vector<LightCone*> {
 public:
  ~LightCones();
  void resize(const size_type n);
  void clear();
  
};

#endif
