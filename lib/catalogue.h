#ifndef CATALOGUE
#define CATALOGUE 1

#include <vector>
#include <deque>
#include "particle.h"
#include "kdtree.h"
#include "hod.h"
#include "lightcone.h"

class Catalogue: public std::vector<Particle>
{
 public:
  Catalogue();
  
  KDTree* tree;
  size_t ntree;
  int ncen, nsat;
};

class Catalogues: public std::deque<Catalogue*> {
 public:
  Catalogues();
  Catalogues(const size_t ncat);
  ~Catalogues();
  void allocate(const size_t ncat);
};

void catalogue_init();
void catalogue_free();

void catalogue_generate_mock(Hod* const hod,
			     LightCone const * const lightcone,
			     const double z_min, const double z_max,
			     Catalogue * const cat);

void catalogue_generate_random(Hod* const hod,
			       LightCone const * const lightcone,
			       const double z_min, const double z_max,
			       Catalogue * const cat);


#endif
