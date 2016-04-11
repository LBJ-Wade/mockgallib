#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "power.h"
#include "msg.h"

PowerSpectrum* power_alloc(const char filename[])
{
  // Read k P from filename
  //PowerSpectrum* const ps= (PowerSpectrum*) malloc(sizeof(PowerSpectrum)); assert(ps);
  PowerSpectrum* const ps= new PowerSpectrum();

  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    msg_printf(msg_error, "Error: Unable to open input power spectrum file: %s\n",filename);
    throw PowerFileError();
  }

  int nalloc= 1000;
  double* buf= (double*) malloc(sizeof(double)*2*nalloc);

  char line[128];
  int nlines= 0;
  double k, P;

  // Read lines and push to buf as k1,P1,k2,P2, ...
  // Doubles the length of buf when the length of the array is not enough
  while(fgets(line, 127, fp)) {
    if(nlines == nalloc) {
      msg_printf(msg_debug, "Reallocating power spectrum table %d -> %d\n",
		 nalloc, 2*nalloc);
      nalloc *= 2;
      buf= (double*) realloc(buf, sizeof(double)*2*nalloc); assert(buf);
    }
    
    if(line[0] == '#')
      continue;
    else if(sscanf(line, "%lg %lg", &k, &P) == 2) {
      buf[2*nlines    ]= k;
      buf[2*nlines + 1]= P;
      
      nlines++;
    }
    else {
      msg_printf(msg_warn,
	    "Warning: Unable to understand a line in the power spectrum file;"
	    "following data are ignored: %s", line);
      break;
    }
  }

  int ret= fclose(fp); assert(ret == 0);
  
  msg_printf(msg_verbose, "Read %d pairs of k P(k) from %s\n", nlines, filename);

  ps->k   = (double*) malloc(2*nlines*sizeof(double)); assert(ps->k);
  ps->P   = ps->k + nlines;
  
  for(int j=0; j<nlines; j++) {
    ps->k[j] = buf[2*j];
    ps->P[j] = buf[2*j + 1];
  }
  free(buf);
  
  ps->n= nlines;
  
  return ps;
}

void power_free(PowerSpectrum* const ps)
{
  free(ps->k);
  delete ps;
}  

double power_sigma(PowerSpectrum const * const ps, const double R)
{
  // Computes sigma (rms amplituede) smoothed on scale R
  // R: smoothing length [/h Mpc] (8 for sigma_8)
  // 1/(2*pi^2) \int P(k) W(k*R) k^2 dk
  
  const double fac= 1.0/(2.0*M_PI*M_PI);

  double k0= ps->k[0];
  double f0= 0.0;
  
  double sigma2= 0.0;
  for(int i=0; i<ps->n; i++) {
    double k= ps->k[i];
    double x= k*R;
    double w= 3.0*(sin(x)-x*cos(x))/(x*x*x);
    double f1= ps->P[i]*k*k*w*w;
    
    sigma2 += 0.5*(f0 + f1)*(k - k0);

    k0= k;
    f0= f1;
  }

  return sqrt(fac*sigma2);
}
