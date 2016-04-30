#include <vector>
#include <cassert>
#include "msg.h"
#include "util.h"
#include "halo_mass.h"

using namespace std;

static void read_nfof_mass(const char filename[],
			   vector<double>& v_nfof, vector<double>& v_M);

HaloMassFoF::HaloMassFoF(const char filename[], const double M0) :
  M0_(M0)
{
  // Reads nfof -> M file
  
  vector<double> v_nfof, v_M;

  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    msg_printf(msg_fatal,
	       "Error: unable to open FoF -> M calibration file %s\n",
	       filename);
    throw ErrorFile();
  }

  char buf[128];
  while(fgets(buf, 127, fp)) {
    if(buf[0] == '#')
      continue;

    int nfof;
    double M;

    int ret= sscanf(buf, "%d %le\n", &nfof, &M); assert(ret == 2);
    v_nfof.push_back((double) nfof);
    v_M.push_back(M);

    if(M > M0_) {
      m_ = M/nfof;
      nfof0_= nfof;
      break;
    }
  }

  fclose(fp);

  spline_ = gsl_spline_alloc(gsl_interp_cspline, v_nfof.size());
  acc_ = gsl_interp_accel_alloc();
  gsl_spline_init(spline_, &v_nfof.front(), &v_M.front(), v_nfof.size());

}

HaloMassFoF::~HaloMassFoF()
{
  gsl_spline_free(spline_);
  gsl_interp_accel_free(acc_);
}

