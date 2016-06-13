#include "hod.h"

// x= z - 0.5
// logMmin = c0 + c1*x + c2*x^2 + c3*x^3
// sigma   = c4 + c5*x
// M0      = Mmin
// logM1   = c6 + c7*x
// alpha   = c8 + c9*x

Hod::Hod()
{
  // Initial guess
  c[0]= 12.0;
  c[1]= 0.0;
  c[2]= 2.0;
  c[3]= 0.0;
  c[4]= 0.1;
  c[5]= 0.0;
  c[6]= 15.0;
  c[7]= 0.0;
  c[8]= 1.5;
  c[9]= 0.0;
}

const double Hod::z0= 0.5;


