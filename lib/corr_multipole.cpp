#ifdef WITHMPI
#include <mpi.h>
#endif

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "msg.h"
#include "comm.h"
#include "catalogue.h"
#include "corr_multipole.h"

using namespace std;

static float rmax2;
static int nbin, nbin_mu;
static float r_min, r_max;
static KDTree* tree_alloc= 0;
static vector<CorrMultipole*> vcorr;
static float ra_min=0.0, dec_min= 0.0;

void corr_multipole_init(const float r_min_, const float r_max_, const int nbin_, const int nbin_mu_)
{
  r_min= r_min_;
  r_max= r_max_;
  nbin= nbin_;
  nbin_mu= nbin_mu_;
}


void corr_multipole_free()
{
  free(tree_alloc);
}


//
// Local function declairations
//
static void allocate_vcorr(const size_t n_data_cat);

static size_t count_pairs_auto(KDTree const * const tree,
			const size_t ntree,
			Histogram2D<LogBin, LinearBin>* const hist);

static size_t count_pairs_cross(KDTree const * const tree1, const size_t ntree1,
			 KDTree const * const tree2,
			 Histogram2D<LogBin, LinearBin>* const hist);

static void compute_corr_from_histogram2d(
			 Histogram2D<LogBin, LinearBin> const * const dd,
			 Histogram2D<LogBin, LinearBin> const * const dr,
			 Histogram2D<LogBin, LinearBin> const * const rr,
			 CorrMultipole* const corr);

static void accumulate_hist(Histogram2D<LogBin, LinearBin>* const hist);
static size_t count_num_points(Catalogues const * const v);

//
// Inline helper functions
//
static inline float dist1(const float left1, const float right1,
			  const float left2, const float right2)
{
  // Distance between two segments [left1,right1] and [left2,right2]
  if(right1 < left2)
    return left2 - right1;
  else if(right2 < left1)
    return left1 - right2;

  return 0.0f;
}


static inline float sq(const float x[])
{
  return (double) x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}


static inline float norm(const float x[])
{
  return sqrt((double) x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}


static inline void dist_spherical(const float x[], const float y[], float& r, float& mu)
{
  // pi = (x - y).\hat{(x + y)/2}
  //    = (|x|^2 - |y|^2)/|x + y|

  float dx[3];
  dx[0]= x[0] - y[0];
  dx[1]= x[1] - y[1];
  dx[2]= x[2] - y[2];
  float r2= dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
  r= sqrt(r2);

  // Line-of-sight unit vector \hat{(x+y)/2}
  float rhat[3];
  rhat[0]= x[0] + y[0];
  rhat[1]= x[1] + y[1];
  rhat[2]= x[2] + y[2];
  float sum_norm= norm(rhat); assert(sum_norm > 0.0f);
  rhat[0] /= sum_norm;
  rhat[1] /= sum_norm;
  rhat[2] /= sum_norm;

  // You can expand (x-y).(x+y) = |x|^2 - |y|^2 but it introduces larger
  // round-off error
  
  float pi= fabs(dx[0]*rhat[0] + dx[1]*rhat[1] + dx[2]*rhat[2]);
  mu= pi/r;
}


size_t count_num_points(Catalogues const * const v)
{
  // Counts total number of points in the catalogues
  size_t n=0;
  
  for(Catalogues::const_iterator cat= v->begin(); cat != v->end(); ++cat) {
    n += (*cat)->size();
  }
  return n;
}

void corr_multipole_radec_min(const float ra_min_, const float dec_min_)
{
  ra_min= ra_min_;
  dec_min= dec_min_;

  msg_printf(msg_verbose, "set corr_multipole ra-dec min %e %e\n", ra_min, dec_min);
}



static size_t count_pairs_leaf_tree_auto(KDTree const * const leaf,
					 KDTree const * const tree,
  				    Histogram2D<LogBin, LinearBin>* const hist)
{
  // prerequisit: set rmax
  //            : binary tree build
  assert(leaf->subtree[0] == 0 && leaf->subtree[1] == 0);
	 
  if(tree - leaf > 0)
    return 0;
  // No need to double count
  // By construction subtree index is always larger than the tree index
  
  float dx= dist1(leaf->left[0], leaf->right[0],
		  tree->left[0], tree->right[0]);
  float dy= dist1(leaf->left[1], leaf->right[1],
		  tree->left[1], tree->right[1]);
  float dz= dist1(leaf->left[2], leaf->right[2],
  		  tree->left[2], tree->right[2]);

  const float r2= dx*dx + dy*dy + dz*dz;

  // leaf and tree are far enough
  if(r2 > rmax2)
    return 0;

  size_t count= 0;    

  if(tree->subtree[0] == 0 && tree->subtree[1] == 0) {
    // This tree is a leaf (no further subtree)
    // leaf - leaf pair count

    float r, mu;

    if(leaf == tree) {
      for(Particle const *p= leaf->particles[0]; p != leaf->particles[1]; ++p) {
	for(Particle const *q= p+1; q != leaf->particles[1]; ++q) {
	  if(fabs(p->radec[0] - q->radec[0]) >= ra_min ||
	     fabs(p->radec[1] - q->radec[1]) >= dec_min){

	    dist_spherical(p->x, q->x, r, mu);
	  
	    count++;
	    hist->add(r, mu, p->w * p->w);
	  }
	}
      }
    }
    else {
      for(Particle const *p= leaf->particles[0]; p != leaf->particles[1]; ++p){
	for(Particle const *q= tree->particles[0]; q != tree->particles[1];++q){
	  if(fabs(p->radec[0] - q->radec[0]) >= ra_min ||
	     fabs(p->radec[1] - q->radec[1]) >= dec_min){

	    dist_spherical(p->x, q->x, r, mu);
	  
	    count++;
	    hist->add(r, mu, p->w * q->w);
	  }
	}
      }
    }

    return count;
  }

  // Recursively seach subtree
  if(tree->subtree[0])
    count += count_pairs_leaf_tree_auto(leaf, tree->subtree[0], hist);
  if(tree->subtree[1])
    count += count_pairs_leaf_tree_auto(leaf, tree->subtree[1], hist);

  return count;
}


size_t count_pairs_auto(KDTree const * const tree,
			const size_t ntree,
			Histogram2D<LogBin, LinearBin>* const hist)
{
  // Run count_pairs_leaf_tree for each leaf
  size_t count= 0;

  for(size_t i=0; i<ntree; ++i) {
    KDTree const * leaf= tree + i;
    if(leaf->subtree[0] == 0 && leaf->subtree[1] == 0) {
      count += count_pairs_leaf_tree_auto(leaf, tree, hist);
    }
  }

  return count;
}


static size_t count_pairs_leaf_tree_cross(KDTree const * const leaf,
				   KDTree const * const tree,
				   Histogram2D<LogBin, LinearBin>* const hist)
{
  // prerequisit: set rmax
  //            : binary tree build
  assert(leaf->subtree[0] == 0 && leaf->subtree[1] == 0);
	 
  float dx= dist1(leaf->left[0], leaf->right[0],
		  tree->left[0], tree->right[0]);
  float dy= dist1(leaf->left[1], leaf->right[1],
		  tree->left[1], tree->right[1]);
  float dz= dist1(leaf->left[2], leaf->right[2],
  		  tree->left[2], tree->right[2]);

  const float r2= dx*dx + dy*dy + dz*dz;

  // leaf and tree are far enough
  if(r2 > rmax2)
    return 0;

  size_t count= 0;    

  if(tree->subtree[0] == 0 && tree->subtree[1] == 0) {
    // This tree is a leaf (no further subtree)
    // leaf - leaf pair count

    float r, mu;

    for(Particle const *p= leaf->particles[0]; p != leaf->particles[1]; ++p){
      for(Particle const *q= tree->particles[0]; q != tree->particles[1];++q){
	if(fabs(p->radec[0] - q->radec[0]) >= ra_min ||
	   fabs(p->radec[1] - q->radec[1]) >= dec_min){
	  dist_spherical(p->x, q->x, r, mu);
	
	  count++;
	  hist->add(r, mu, p->w * q->w);
	}
      }
    }

    return count;
  }

  // Recursively seach subtree
  if(tree->subtree[0])
    count += count_pairs_leaf_tree_cross(leaf, tree->subtree[0], hist);
  if(tree->subtree[1])
    count += count_pairs_leaf_tree_cross(leaf, tree->subtree[1], hist);

  return count;
}


size_t count_pairs_cross(KDTree const * const tree1, const size_t ntree1,
			 KDTree const * const tree2,
			 Histogram2D<LogBin, LinearBin>* hist)
{
  // Run count_pairs_leaf_tree for each leaf
  size_t count= 0;

  for(size_t i=0; i<ntree1; ++i) {
    KDTree const * leaf= tree1 + i;
    if(leaf->subtree[0] == 0 && leaf->subtree[1] == 0) {
      count += count_pairs_leaf_tree_cross(leaf, tree2, hist);
    }
  }

  return count;
}


//
// Struct Multipole
//
CorrMultipole::CorrMultipole(const int nbin) :
  n(nbin)
{
  r= (double*) malloc(sizeof(double)*nbin*3);
  xi0= r + nbin;
  xi2= r + nbin;
}

CorrMultipole::~CorrMultipole()
{
  free(r);
}

void allocate_vcorr(const size_t n_data_cat)
{
  assert(nbin > 0);
  if(vcorr.size() == n_data_cat)
    return;


  for(vector<CorrMultipole*>::iterator
	corr= vcorr.begin(); corr != vcorr.end(); ++corr)
    delete *corr;
  vcorr.clear();
  
  for(int i=0; i<n_data_cat; ++i) {
    CorrMultipole* const corr= new CorrMultipole(nbin);
    vcorr.push_back(corr);
  }
}

CorrMultipole* corr_multipole_i(const int i)
{
  return vcorr.at(i);
}


void compute_corr_from_histogram2d(
	   Histogram2D<LogBin, LinearBin> const * const dd,
	   Histogram2D<LogBin, LinearBin> const * const dr,
	   Histogram2D<LogBin, LinearBin> const * const rr,
	   CorrMultipole* const corr)
{
  // Integrate 2D historgram xi(r, mu) = (DD - 2 DR + RR)/RR to multipoles

  const int nx= dd->x_nbin();
  const int ny= dd->y_nbin();
  const double dmu= 1.0/ny;
  
  assert(rr->npairs > 0);
  assert(corr->n == nx);

  for(int ix=0; ix<nx; ++ix) {
    corr->r[ix]= dd->x_bin(ix);
    double xi0= 0.0;
    double xi2= 0.0;
    
    for(int iy=0; iy<ny; ++iy) {
      int index= ix*dd->y_nbin() + iy;
      double mu = dd->y_bin(iy);

      double rrr= rr->hist[index]/rr->npairs; assert(rr->npairs > 0);
      if(rrr > 0.0) {
	double xi = (dd->hist[index]/dd->npairs - 2.0*dr->hist[index]/dr->npairs)/rrr + 1.0;
	xi0 += xi*dmu;
	xi2 += (4.5*mu*mu - 1.5)*xi*dmu;
      }
      // xi = (DD - 2*DR + RR)/RR
      // xi[l] = (2l + 1) int_0^1 P_l(mu) xi(r, mu) dmu
    }    
    corr->xi0[ix]= xi0;
    corr->xi2[ix]= xi2;
	  
    assert(!isnan(xi0));
    assert(!isnan(xi2));
  }
}


void accumulate_hist(Histogram2D<LogBin, LinearBin>* const hist)
{
  if(comm_n_nodes() == 1)
    return;

  //
  // Accumulate from all MPI nodes
  //
  
  const int nx= hist->x_nbin();
  const int ny= hist->y_nbin();
  const int n=  hist->size();
  
  double npairs= hist->npairs;
  
#ifdef WITHMPI
  double* const hist_sum= (double*) malloc(sizeof(double)*n); assert(hist_sum);
  double npairs_sum;

  MPI_Allreduce(&npairs, &npairs_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  msg_printf(msg_debug, "Accumulate hist %.1lf pairs -> %lf; %.3f\n",
	     npairs, npairs_sum, npairs_sum / npairs);

  MPI_Allreduce(hist->hist, hist_sum, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  hist->npairs= npairs_sum;
  for(int i=0; i<n; ++i){
    hist->hist[i]= hist_sum[i];
  }

  free(hist_sum);
#endif
}

		     
void corr_multipole_compute_pairs_rr(Catalogues* const cats_rand,
				   Histogram2D<LogBin, LinearBin>* const rr)
{
  rmax2= r_max*r_max;

  //
  // Setup KDTree
  //
  //const int n= cats->size();
  const int quota = 32;
    
  if(tree_alloc == 0) {
    size_t nalloc= count_num_points(cats_rand);


    tree_alloc= (KDTree*) malloc(sizeof(KDTree)*nalloc);
    msg_printf(msg_verbose, "%lu trees allocated (%lu Mbytes)\n",
	       nalloc, nalloc*sizeof(KDTree) / (1024*1024));
  }

  KDTree* tree_free= tree_alloc;
  size_t ntree_used= 0;

  // KDTree for Randoms
  for(Catalogues::iterator cat= cats_rand->begin();
      cat != cats_rand->end(); ++cat) {
    (*cat)->tree= tree_free;
    (*cat)->ntree = kdtree_construct((*cat)->tree, &((*cat)->front()),
				     (*cat)->size(), quota);
    ntree_used += (*cat)->ntree;
    tree_free += (*cat)->ntree;
  }
  
  msg_printf(msg_verbose, "%lu trees used (%lu Mbytes).\n",
	     ntree_used, ntree_used*sizeof(KDTree) / (1024*1024));
  msg_printf(msg_verbose, "Count RR pairs.\n");

  //
  // Count pairs RR
  //
  rr->clear();
  for(Catalogues::iterator cat=
	cats_rand->begin(); cat != cats_rand->end(); ++cat) {    
    if(!(*cat)->empty()) {
      count_pairs_auto((*cat)->tree, (*cat)->ntree, rr);  
      rr->npairs += 0.5*(*cat)->size()*((*cat)->size()-1);
    }
  }

  msg_printf(msg_verbose, "%.1lf RR pairs.\n", rr->npairs);

  accumulate_hist(rr);
}

void corr_multipole_compute_with_rr(Catalogues* const cats_data,
				    Catalogues* const cats_rand,
				    Histogram2D<LogBin, LinearBin> const * const rr)
{
  // Setup vcorr
  allocate_vcorr(cats_data->size());

  rmax2= r_max*r_max;

  assert(nbin == rr->x_nbin());
  assert(nbin_mu == rr->y_nbin());

  //
  // Setup KDTree
  //
  const int quota = 32;
    
  if(tree_alloc == 0) {
    size_t nalloc= count_num_points(cats_data) + count_num_points(cats_rand);

    tree_alloc= (KDTree*) malloc(sizeof(KDTree)*nalloc);
    msg_printf(msg_verbose, "%lu trees allocated (%lu Mbytes)\n",
	       nalloc, nalloc*sizeof(KDTree) / (1024*1024));
  }

  KDTree* tree_free= tree_alloc;
  size_t ntree_used= 0;

  // KDTree for data
  for(Catalogues::iterator cat= cats_data->begin();
      cat != cats_data->end(); ++cat) {
    (*cat)->tree= tree_free;
    (*cat)->ntree = kdtree_construct((*cat)->tree, &((*cat)->front()),
				     (*cat)->size(), quota);
    ntree_used += (*cat)->ntree;
    tree_free += (*cat)->ntree;
  }

  // KDTree for Randoms
  for(Catalogues::iterator cat= cats_rand->begin();
      cat != cats_rand->end(); ++cat) {
    (*cat)->tree= tree_free;
    (*cat)->ntree = kdtree_construct((*cat)->tree, &((*cat)->front()),
				     (*cat)->size(), quota);
    ntree_used += (*cat)->ntree;
    tree_free += (*cat)->ntree;
  }
  
  msg_printf(msg_verbose, "%lu trees used (%lu Mbytes).\n",
	     ntree_used, ntree_used*sizeof(KDTree) / (1024*1024));

  msg_printf(msg_verbose, "Count DD and DR pairs.\n");

  Histogram2D<LogBin, LinearBin>
    dd(LogBin(r_min, r_max, nbin), LinearBin(0.0f, 1.0f, nbin_mu)),
    dr(LogBin(r_min, r_max, nbin), LinearBin(0.0f, 1.0f, nbin_mu));

  int icat= 0;
  for(Catalogues::iterator cat= cats_data->begin();
      cat != cats_data->end(); ++cat) {
    dd.clear();
    dr.clear();

    // DD
    if(!(*cat)->empty())
      count_pairs_auto((*cat)->tree, (*cat)->ntree, &dd);

    dd.npairs += 0.5*(*cat)->size()*((*cat)->size()-1);


    // DR for all random catalogues
    for(Catalogues::iterator rcat= cats_rand->begin();
      rcat != cats_rand->end(); ++rcat) {    

      if(!(*cat)->empty() && !(*rcat)->empty())
	count_pairs_cross((*cat)->tree, (*cat)->ntree, (*rcat)->tree, &dr);
      dr.npairs += (*cat)->size()*(*rcat)->size();
    }

    accumulate_hist(&dd);
    accumulate_hist(&dr);

    // compute correlation function
    compute_corr_from_histogram2d(&dd, &dr, rr, vcorr.at(icat));

    icat++;
  }

}
