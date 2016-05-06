#include <cassert>
#include "lightcone.h"

using namespace std;

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

void LightCones::resize(const size_type n)
{
  for(size_type i=size(); i<n; ++i)
    push_back(new LightCone());

  assert(size() >= n);
}
