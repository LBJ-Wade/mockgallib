#include "lightcone.h"

LightCones::~LightCones()
{
  for(LightCones::iterator p= begin(); p != end(); ++p) {
    delete *p;
  }
}
