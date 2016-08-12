#ifndef HIST2D_H
#define HIST2D_H 1

#include <cstdlib>

struct LinearBin {
  LinearBin(const float x_min, const float x_max, const int n_bin) :
    xmin(x_min), xmax(x_max), nbin(n_bin)
  {
    dx_inv= nbin/(xmax - xmin);
  }
  int operator()(const float x) const {
    //if(xmin < x && x <= xmax) td::cerr << "y OK\n";
    // return the index of the bin -- can be outside the range [0,nbin)
    return (int)floor((x - xmin)*dx_inv);
  }
  float x_bin(const int ix, const float offset=0.5f) const {
    return xmin + (xmax - xmin)/nbin*(ix + offset);
  }
  float xmin, xmax, dx_inv;
  int nbin;
};

struct LogBin {
LogBin(const float x_min, const float x_max, const int n_bin) :
  nbin(n_bin)
  {
    xmin= log(x_min);
    xmax= log(x_max);
    dx_inv= nbin/(xmax - xmin);
  }
  int operator()(const float x) const {
    //if(xmin < log(x) && log(x) < xmax) std::cerr << "x OK";
    //fprintf(stderr, "logbin %e %e %e %e\n", x, log(x), xmin, xmax);
    return (int)floor((log(x) - xmin)*dx_inv);
  }
  float x_bin(const int ix, const float offset=0.5f) const {
    return exp(xmin + (xmax - xmin)/nbin*(ix + offset));
  }

  float xmin, xmax, dx_inv;
  int nbin;
};

template<class X, class Y>
class Histogram2D {
 public:
  Histogram2D(X xbin_, Y ybin_);
  ~Histogram2D();
  void add(float x, float y, float val) {
    int ix= xbin(x);
    int iy= ybin(y);
    if(0 <= ix && ix < xbin.nbin && 0 <= iy && iy < ybin.nbin) {
      int index= ix*ybin.nbin + iy;
#pragma omp atomic
      hist[index] += val;
    }
  }
  int size() const {
    return xbin.nbin * ybin.nbin;
  }
  int x_nbin() const {
    return xbin.nbin;
  }
  int y_nbin() const {
    return ybin.nbin;
  }
  double total() const {
    double sum= 0.0;
    const int n= xbin.nbin*ybin.nbin;
    for(int i=0; i<n; ++i) sum += hist[i];
    return sum;
  }
  double operator[](const int i) const {
    return hist[i];
  }
  float x_bin(const int ix, const float offset=0.5f) const {
    return xbin.x_bin(ix, offset);
  }
  float y_bin(const int iy, const float offset=0.5f) const {
    return ybin.x_bin(iy, offset);
  }
  double operator()(const int ix, const int iy) const {
    return hist[ix*ybin.nbin + iy];
  }

  void operator+=(const Histogram2D<X,Y>& h) {
    const int n= size(); assert(n == h.size());
    for(int i=0; i<n; ++i)
      hist[i] += h.hist[i];
  }
  void operator-=(const Histogram2D<X,Y>& h) {
    const int n= size(); assert(n == h.size());
    for(int i=0; i<n; ++i)
      hist[i] -= h.hist[i];
  }
  void clear() {
    const int n= size();
    for(int i=0; i<n; ++i)
      hist[i] = 0.0;
  }

  double* hist;
 private:
  X xbin;
  Y ybin;
};

template<class X, class Y>
  Histogram2D<X,Y>::Histogram2D(X xbin_, Y ybin_) : xbin(xbin_), ybin(ybin_)
{
  //const int n= xbin.nbin*ybin.nbin];
  //hist= new double[xbin.nbin*ybin.nbin];
  hist= (double*) calloc(sizeof(double), xbin.nbin*ybin.nbin);
}

template<class X, class Y>
  Histogram2D<X,Y>::~Histogram2D()
{
  free(hist);
  //delete [] hist;
}


#endif
