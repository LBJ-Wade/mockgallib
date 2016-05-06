#include "lightcone.h"

LightCones::~LightCones()
{
  for(LightCones::iterator p= begin(); p != end(); ++p) {
    delete *p;
  }
}

void LightCones::clear()
{
  for(LightCones::iterator p= begin(); p != end(); ++p) {
    (*p)->clear();
  }
}
