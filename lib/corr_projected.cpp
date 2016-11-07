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
#include "corr_projected.h"
#include "hist2d.h"
#include "corr_pair_correction.h"

using namespace std;

static float rmax2;
static int nbin, nbin_pi;
static float rp_min, rp_max, pi_max;
static KDTree* tree_alloc= 0;
static vector<CorrProjected*> vcorr;
static float ra_min=0.0, dec_min= 0.0;
static bool pair_correction= false;

//define DEBUG_TREE 1  // disable kdtree and count N^2 pairs 

void corr_projected_init(const float rp_min_, const float rp_max_, const int nbin_, const float pi_max_, const int nbin_pi_)
{
  rp_min= rp_min_;
  rp_max= rp_max_;
  nbin=   nbin_;
  pi_max= pi_max_;
  nbin_pi= nbin_pi_;

  rmax2= rp_max*rp_max + pi_max*pi_max;

  msg_printf(msg_debug,
	     "corr_projected_init %e %e %d %e %d\n",
	     rp_min, rp_max, nbin, pi_max, nbin_pi);
}

void corr_set_pair_correction(const char filename[])
{
  corr_pair_correction_init(filename);
  pair_correction= true;
}


void corr_projected_free()
{
  free(tree_alloc);
}


//
// Local function declairations
//
static void allocate_vcorr(const size_t n_data_cat);

static size_t count_pairs_auto(KDTree const * const tree,
			       const size_t ntree,
			       const bool is_dd,
			Histogram2D<LogBin, LinearBin>* const hist);

static size_t count_pairs_cross(KDTree const * const tree1, const size_t ntree1,
			 KDTree const * const tree2,
			 Histogram2D<LogBin, LinearBin>* const hist);

static void count_pairs_auto_direct(Catalogue const * const cat,
			     const bool is_dd,
			     Histogram2D<LogBin, LinearBin>* const hist);

static void count_pairs_cross_direct(Catalogue const * const cat,
			      Catalogue const * const rcat,
			      Histogram2D<LogBin, LinearBin>* const hist);

static void compute_corr_from_histogram2d(
			 Histogram2D<LogBin, LinearBin> const * const dd,
			 Histogram2D<LogBin, LinearBin> const * const dr,
			 Histogram2D<LogBin, LinearBin> const * const rr,
			 const double pi_max,
			 CorrProjected* const corr);

static void accumulate_hist(Histogram2D<LogBin, LinearBin>* const hist);
static void corr_projected_summarise(CorrProjected* const corr);
static size_t count_num_points(Catalogues const * const v);

void compute_wsum(Catalogue * const cat);

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


static inline void dist_cylinder(const float x[], const float y[], float& rp, float& pi)
{
  // pi = (x - y).\hat{(x + y)/2}
  //    = (|x|^2 - |y|^2)/|x + y|

  float dx[3];
  dx[0]= x[0] - y[0];
  dx[1]= x[1] - y[1];
  dx[2]= x[2] - y[2];
  float r2= dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

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
  
  pi= fabs(dx[0]*rhat[0] + dx[1]*rhat[1] + dx[2]*rhat[2]);

  if(r2 < pi*pi) {
    float err=(r2 - pi*pi) / r2;
    assert(fabs(err) < 1.0e-5);
    rp= 0.0f;
  }
  else
    rp= sqrt(r2 - pi*pi);
}

static inline float dist_angle(const float radec1[], const float radec2[])
{
  const float dra= radec1[0] - radec2[0];
  const float ddec= radec1[1] - radec2[1];
  
  return sqrt(dra*dra + ddec*ddec);
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

void set_radec_min(const float ra_min_, const float dec_min_)
{
  ra_min= ra_min_;
  dec_min= dec_min_;

  msg_printf(msg_verbose, "ra-dec min %e %e\n", ra_min, dec_min);
}

    
void corr_projected_compute(Catalogues* const cats_data,
			    Catalogues* const cats_rand,
			    CorrProjected* const corr)
{
  cerr << "waiting at a barrier corr_projected_compute\n";
  comm_barrier();
  //
  // Compute 2D correlation function
  //
  rmax2= rp_max*rp_max + pi_max*pi_max;

  // Setup vcorr
  allocate_vcorr(cats_data->size());
  
  //
  // Setup KDTree
  //
  size_t nalloc= count_num_points(cats_data) + count_num_points(cats_rand);
  const int quota = 32;

  if(tree_alloc == 0) {
    tree_alloc= (KDTree*) malloc(sizeof(KDTree)*nalloc);
    msg_printf(msg_verbose, "%lu trees allocated (%lu Mbytes)\n",
	       nalloc, nalloc*sizeof(KDTree) / (1024*1024));
  }

  KDTree* tree_free= tree_alloc;
  size_t ntree_used= 0;

  // KDTree for Data
  for(Catalogues::iterator cat= cats_data->begin();
      cat != cats_data->end(); ++cat) {
    (*cat)->tree= tree_free;

    (*cat)->ntree = kdtree_construct((*cat)->tree, &((*cat)->front()),
				     (*cat)->size(), quota);
    ntree_used += (*cat)->ntree;
    tree_free +=  (*cat)->ntree;
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
  msg_printf(msg_verbose, "Computing correlation function.\n");



  //
  // Setup Histogram
  //
  Histogram2D<LogBin, LinearBin>
    dd(LogBin(rp_min, rp_max, nbin), LinearBin(0.0f, pi_max, nbin_pi)),
    dr(LogBin(rp_min, rp_max, nbin), LinearBin(0.0f, pi_max, nbin_pi)),
    rr(LogBin(rp_min, rp_max, nbin), LinearBin(0.0f, pi_max, nbin_pi));

  //
  // Count pairs
  //

  // RR
  size_t debug_count_dd= 0, debug_count_dr=0, debug_count_rr=0;
  size_t count_dd=0, count_dr=0, count_rr=0;
  
  rr.clear();
  for(Catalogues::iterator rcat= cats_rand->begin();
      rcat != cats_rand->end(); ++rcat) {    
    if(!(*rcat)->empty()) {
      compute_wsum(*rcat);
      debug_count_rr +=
	count_pairs_auto((*rcat)->tree, (*rcat)->ntree, false, &rr);
      rr.npairs += 0.5*((*rcat)->wsum*(*rcat)->wsum - (*rcat)->w2sum);
    }
    
    count_rr += (*rcat)->size()*((*rcat)->size() - 1)/2;
  }

  accumulate_hist(&rr);

  int icat=0;
  assert(cats_data->size() > 0);
  const int rand_cats_factor= cats_rand->size() / cats_data->size();
  
  for(Catalogues::iterator cat= cats_data->begin();
      cat != cats_data->end(); ++cat) {
    dd.clear();
    dr.clear();

    // DD
    if(!(*cat)->empty()) {
      debug_count_dd +=
	count_pairs_auto((*cat)->tree, (*cat)->ntree, true, &dd);

      compute_wsum(*cat);
      dd.npairs += 0.5*((*cat)->wsum*(*cat)->wsum - (*cat)->w2sum);
    }

    count_dd += (*cat)->size()*((*cat)->size()-1)/2;


    // DR
    // Not ndata_cat * nrand_cat crosses, but nrand_cat cross pairs
    for(size_t icat_rand=0; icat_rand<rand_cats_factor; ++icat_rand) {
      Catalogues::iterator rcat= cats_rand->begin() + icat +
	                           icat_rand * cats_data->size();
      assert(icat + icat_rand * cats_data->size() < cats_rand->size());

      if(!(*cat)->empty() && !(*rcat)->empty()) {
	debug_count_dr += count_pairs_cross((*cat)->tree, (*cat)->ntree, (*rcat)->tree, &dr);
      
	dr.npairs += (*cat)->wsum * (*rcat)->wsum;
	count_dr += (*cat)->size()*(*rcat)->size();
      }
    }

#ifdef DEBUG_TREE
    msg_printf(msg_debug, "debug_count dd %llu %llu\n",
	       debug_count_dd, count_dd);
    msg_printf(msg_debug, "debug_count dr %llu %llu\n",
	       debug_count_dr, count_dr);
    msg_printf(msg_debug, "debug_count rr %llu %llu\n",
	       debug_count_rr, count_rr);

    assert(debug_count_dd == count_dd);
    assert(debug_count_dr == count_dr);
    assert(debug_count_rr == count_rr);
#endif

    accumulate_hist(&dd);
    accumulate_hist(&dr);

    /*
    // ndata_cat * nrand_cat cross pairs
    for(Catalogues::iterator rcat= cats_rand->begin();
      rcat != cats_rand->end(); ++rcat) {    

      if(!(*cat)->empty() && !(*rcat)->empty())
	count_dr_debug += count_pairs_cross((*cat)->tree, (*cat)->ntree, (*rcat)->tree, dr);
      npairs_DR += (*cat)->size()*(*rcat)->size();
    }
    */

    compute_corr_from_histogram2d(&dd, &dr, &rr, pi_max, vcorr.at(icat));

    icat++;
  }
  corr_projected_summarise(corr);
}


static size_t count_pairs_leaf_tree_auto(KDTree const * const leaf,
					 KDTree const * const tree,
					 const bool is_dd,
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
#ifndef DEBUG_TREE
  if(r2 > rmax2)
    return 0;
#endif

  size_t count= 0;
  double pw= 1.0; // pair-wise weight

  if(tree->subtree[0] == 0 && tree->subtree[1] == 0) {
    // This tree is a leaf (no further subtree)
    // leaf - leaf pair count

    float rp, pi;

    if(leaf == tree) {
      for(Particle const *p= leaf->particles[0]; p != leaf->particles[1]; ++p) {
	for(Particle const *q= p+1; q != leaf->particles[1]; ++q) {
	  count++;
	  if(fabs(p->radec[0] - q->radec[0]) >= ra_min ||
	     fabs(p->radec[1] - q->radec[1]) >= dec_min) {
	    
	    dist_cylinder(p->x, q->x, rp, pi);

	    if(is_dd && pair_correction)
	      pw= corr_pair_correction(dist_angle(p->radec, q->radec));

	    hist->add(rp, pi, p->w * q->w / pw);
	  }
	}
      }
    }
    else {
      for(Particle const *p= leaf->particles[0]; p != leaf->particles[1]; ++p){
	for(Particle const *q= tree->particles[0]; q != tree->particles[1];++q){
	  count++;
		    
	  if(fabs(p->radec[0] - q->radec[0]) >= ra_min ||
	     fabs(p->radec[1] - q->radec[1]) >= dec_min) {

	    dist_cylinder(p->x, q->x, rp, pi);

	    if(is_dd && pair_correction)
	      pw= corr_pair_correction(dist_angle(p->radec, q->radec));

	    hist->add(rp, pi, p->w * q->w / pw);
	  }
	}
      }
    }

    return count;
  }

  // Recursively seach subtree
  if(tree->subtree[0])
    count += count_pairs_leaf_tree_auto(leaf, tree->subtree[0], is_dd, hist);
  if(tree->subtree[1])
    count += count_pairs_leaf_tree_auto(leaf, tree->subtree[1], is_dd, hist);

  return count;
}


size_t count_pairs_auto(KDTree const * const tree,
			const size_t ntree,
			const bool is_dd,
			Histogram2D<LogBin, LinearBin>* const hist)
{
  // Run count_pairs_leaf_tree for each leaf
  size_t count= 0;

  for(size_t i=0; i<ntree; ++i) {
    KDTree const * leaf= tree + i;
    if(leaf->subtree[0] == 0 && leaf->subtree[1] == 0) {
      count += count_pairs_leaf_tree_auto(leaf, tree, is_dd, hist);
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
#ifndef DEBUG_TREE
  if(r2 > rmax2)
    return 0;
#endif

  size_t count= 0;    

  if(tree->subtree[0] == 0 && tree->subtree[1] == 0) {
    // This tree is a leaf (no further subtree)
    // leaf - leaf pair count

    float rp, pi;

    for(Particle const *p= leaf->particles[0]; p != leaf->particles[1]; ++p){
      for(Particle const *q= tree->particles[0]; q != tree->particles[1];++q){
	count++;

	if(fabs(p->radec[0] - q->radec[0]) >= ra_min ||
	   fabs(p->radec[1] - q->radec[1]) >= dec_min){
	  dist_cylinder(p->x, q->x, rp, pi);

	  hist->add(rp, pi, p->w * q->w);
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
// Struct CorrProjected
//
CorrProjected::CorrProjected(const int nbin) :
  n(nbin)
{
  rp= (double*) malloc(sizeof(double)*nbin*3);
  wp= rp + nbin;
  dwp= wp + nbin;
}

CorrProjected::~CorrProjected()
{
  free(rp);
}

void CorrProjected::print(FILE* fp)
{
  for(int i=0; i<n; ++i) {
    fprintf(fp, "%e %e %e\n", rp[i], wp[i], dwp[i]);
  }
}

void allocate_vcorr(const size_t n_data_cat)
{
  assert(nbin > 0);
  if(vcorr.size() == n_data_cat)
    return;


  for(vector<CorrProjected*>::iterator
	corr= vcorr.begin(); corr != vcorr.end(); ++corr)
    delete *corr;
  vcorr.clear();
  
  for(int i=0; i<n_data_cat; ++i) {
    CorrProjected* const corr= new CorrProjected(nbin);
    vcorr.push_back(corr);
  }
}

CorrProjected* corr_projected_i(const int i)
{
  return vcorr.at(i);
}


void compute_corr_from_histogram2d(
	   Histogram2D<LogBin, LinearBin> const * const dd,
	   Histogram2D<LogBin, LinearBin> const * const dr,
	   Histogram2D<LogBin, LinearBin> const * const rr,
	   const double pi_max,
	   CorrProjected* const corr)
{
  // Project 2D historgram DD, DR, DD
  // to projected correlation function to wp(rp)
  const int nx= dd->x_nbin();
  const int ny= dd->y_nbin();
  const double dpi= 2.0*pi_max / ny; assert(ny > 0);

  assert(rr->npairs > 0);
  assert(corr->n == nx);

  msg_printf(msg_debug,
	     "npairs dd dr rr %e %e %e\n", dd->npairs, dr->npairs, rr->npairs);

  for(int ix=0; ix<nx; ++ix) {
    corr->rp[ix]= dd->x_bin(ix);
    double wp= 0.0;
    for(int iy=0; iy<ny; ++iy) {
      int index= ix*dd->y_nbin() + iy;

      double rrr= rr->hist[index]/rr->npairs; assert(rr->npairs > 0);
      if(rrr > 0.0) {
	wp += ((dd->hist[index]/dd->npairs
	     - 2.0*dr->hist[index]/dr->npairs)/rrr + 1.0)*dpi;
      }

      // xi = (DD - 2*DR + RR)/RR
      // wp = \int_-pi-max^pi-max xi(rp, pi) dpi
    }    
    corr->wp[ix]= wp;
	  
    assert(!isnan(wp));
  }
}


void corr_projected_write(const int index, const vector<CorrProjected*>& vcorr)
{
  assert(!vcorr.empty());
  const int nbin= vcorr.front()->n;
  const int ndat= vcorr.size(); assert(ndat > 0);

  char filename[128];
  sprintf(filename, "wp_%05d.txt", index);
  FILE* fp= fopen(filename, "w"); assert(fp);
  
  for(int i=0; i<nbin; ++i) {
    double rp= vcorr.front()->rp[i];
    double wp_sum= 0.0;
    double wp2_sum= 0.0;
    for(vector<CorrProjected*>::const_iterator corr= vcorr.begin();
	corr != vcorr.end(); ++corr) {
      wp_sum +=  (*corr)->wp[i];
      wp2_sum += (*corr)->wp[i]*(*corr)->wp[i];
    }

    double wp= wp_sum / ndat; assert(ndat > 0);
    double dwp= ndat > 1 ? sqrt((wp2_sum - ndat*wp*wp) / (ndat-1)) : 0.0;

    fprintf(fp, "%e %e %e\n", rp, wp, dwp);
  }
  fclose(fp);
}


void corr_projected_summarise(CorrProjected* const corr)
{
  //
  // Take average of multiple projected correlation functions vcorr
  //
  assert(!vcorr.empty());
  msg_printf(msg_debug, "summarise %lu correlation functions\n",
	     vcorr.size());
  
  const int nbin= vcorr.front()->n;
  const int ndat= vcorr.size(); assert(ndat > 0);

  for(int i=0; i<nbin; ++i) {
    double rp= vcorr.front()->rp[i];
    double wp_sum= 0.0;
    double wp2_sum= 0.0;
    for(vector<CorrProjected*>::const_iterator cp= vcorr.begin();
	cp != vcorr.end(); ++cp) {
      wp_sum +=  (*cp)->wp[i];
      wp2_sum += (*cp)->wp[i]*(*cp)->wp[i];
    }

    double wp= wp_sum / ndat;
    double dwp= ndat > 1 ? sqrt((wp2_sum - ndat*wp*wp)/(ndat-1)) : wp;

    corr->rp[i]= rp;
    corr->wp[i]= wp;
    corr->dwp[i]= dwp;
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

  //cerr << "all reduce accumulate hist\n";
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

		     
void corr_projected_compute_pairs_rr(Catalogues* const cats_rand,
				   Histogram2D<LogBin, LinearBin>* const rr)
{
  msg_printf(msg_verbose, "corr_projected_compute_pairs_rr\n");
  msg_printf(msg_debug, "rp_min= %e, pi_max= %e\n", rp_min, rp_max);

  rmax2= rp_max*rp_max + pi_max*pi_max;

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
      count_pairs_auto((*cat)->tree, (*cat)->ntree, false, rr);

      compute_wsum(*cat);
      rr->npairs += 0.5*((*cat)->wsum*(*cat)->wsum - (*cat)->w2sum);
    }
  }

  msg_printf(msg_verbose, "%.1lf RR pairs.\n", rr->npairs);

  accumulate_hist(rr);
}

void corr_projected_compute_pairs_all(Catalogues* const cats_data,
				     Catalogues* const cats_rand,
				     Histogram2D<LogBin, LinearBin>* const dd,
				     Histogram2D<LogBin, LinearBin>* const dr,
				     Histogram2D<LogBin, LinearBin>* const rr)
{
  msg_printf(msg_verbose, "corr_projected_compute_pairs_all\n");
  msg_printf(msg_debug, "rp_min= %e, pi_max= %e\n", rp_min, rp_max);

  rmax2= rp_max*rp_max + pi_max*pi_max;

  //
  // Setup KDTree
  //
  size_t nalloc= count_num_points(cats_data) + count_num_points(cats_rand);
  const int quota = 32;

  if(tree_alloc == 0) {
    tree_alloc= (KDTree*) malloc(sizeof(KDTree)*nalloc);
    msg_printf(msg_verbose, "%lu trees allocated (%lu Mbytes)\n",
	       nalloc, nalloc*sizeof(KDTree) / (1024*1024));
  }

  KDTree* tree_free= tree_alloc;
  size_t ntree_used= 0;

  // KDTree for Data
  for(Catalogues::iterator cat= cats_data->begin();
      cat != cats_data->end(); ++cat) {
    (*cat)->tree= tree_free;

    (*cat)->ntree = kdtree_construct((*cat)->tree, &((*cat)->front()),
				     (*cat)->size(), quota);
    ntree_used += (*cat)->ntree;
    tree_free +=  (*cat)->ntree;
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
  msg_printf(msg_verbose, "Computing correlation function.\n");

  //
  // Count pairs RR
  //
  rr->clear();
  for(Catalogues::iterator cat=
	cats_rand->begin(); cat != cats_rand->end(); ++cat) {    
    if(!(*cat)->empty()) {
      count_pairs_auto((*cat)->tree, (*cat)->ntree, false, rr);

      compute_wsum(*cat);
      rr->npairs += 0.5*((*cat)->wsum*(*cat)->wsum - (*cat)->w2sum);
    }
  }

  dd->clear();
  dr->clear();

  for(Catalogues::iterator cat= cats_data->begin();
      cat != cats_data->end(); ++cat) {
    // DD
    if(!(*cat)->empty()) {
      count_pairs_auto((*cat)->tree, (*cat)->ntree, true, dd);

      compute_wsum(*cat);
      dd->npairs += 0.5*((*cat)->wsum*(*cat)->wsum - (*cat)->w2sum);
    }

    // DR
    // ndata_cat * nrand_cat cross pairs
    for(Catalogues::iterator rcat= cats_rand->begin();
      rcat != cats_rand->end(); ++rcat) {    

      if(!(*cat)->empty() && !(*rcat)->empty())
	count_pairs_cross((*cat)->tree, (*cat)->ntree, (*rcat)->tree, dr);
      dr->npairs += (*cat)->wsum * (*rcat)->wsum;
    }
  }

  accumulate_hist(dd);
  accumulate_hist(dr);
  accumulate_hist(rr);
}

void corr_projected_compute_with_rr(Catalogues* const cats_data,
				    Catalogues* const cats_rand,
			    Histogram2D<LogBin, LinearBin> const * const rr)
{
  // Setup vcorr
  allocate_vcorr(cats_data->size());

  rmax2= rp_max*rp_max + pi_max*pi_max;

  assert(nbin == rr->x_nbin());
  assert(nbin_pi == rr->y_nbin());

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
    dd(LogBin(rp_min, rp_max, nbin), LinearBin(0.0f, pi_max, nbin_pi)),
    dr(LogBin(rp_min, rp_max, nbin), LinearBin(0.0f, pi_max, nbin_pi));

  int icat= 0;
  for(Catalogues::iterator cat= cats_data->begin();
      cat != cats_data->end(); ++cat) {
    dd.clear();
    dr.clear();

    // DD
    if(!(*cat)->empty())
      count_pairs_auto((*cat)->tree, (*cat)->ntree, true, &dd);

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

    // compute projected correlation function
    compute_corr_from_histogram2d(&dd, &dr, rr, pi_max, vcorr.at(icat));

    icat++;
  }

}

void compute_wsum(Catalogue* const cat)
{
  // cat->wsum  = \sum_i w_i
  // cat->s2sum = \sum_i w_i^2
  double wsum= 0.0, w2sum= 0.0;
  
  for(Catalogue::const_iterator p= cat->begin(); p != cat->end(); ++p) {
    wsum += p->w;
    w2sum += p->w * p->w;
  }

  cat->wsum= wsum;
  cat->w2sum= w2sum;
}


//
// Direct N^2 computation (for tests)
//
void count_pairs_auto_direct(Catalogue const * const cat,
			     const bool is_dd,
			     Histogram2D<LogBin, LinearBin>* const hist)
{
  for(Catalogue::const_iterator p= cat->begin(); p != cat->end(); ++p) {
    for(Catalogue::const_iterator q= p+1; q != cat->end(); ++q) {
      if(fabs(p->radec[0] - q->radec[0]) >= ra_min ||
	 fabs(p->radec[1] - q->radec[1]) >= dec_min) {
	float rp, pi;
	double pw= 1.0;

	dist_cylinder(p->x, q->x, rp, pi);

	if(is_dd && pair_correction)
	  pw= corr_pair_correction(dist_angle(p->radec, q->radec));

	hist->add(rp, pi, p->w * q->w / pw);
      }
    }
  }
}

void count_pairs_cross_direct(Catalogue const * const cat,
			     Catalogue const * const rcat,
			     Histogram2D<LogBin, LinearBin>* const hist)
{
  for(Catalogue::const_iterator p= cat->begin(); p != cat->end(); ++p) {
    for(Catalogue::const_iterator q= rcat->begin(); q != rcat->end(); ++q) {
      if(fabs(p->radec[0] - q->radec[0]) >= ra_min ||
	 fabs(p->radec[1] - q->radec[1]) >= dec_min) {
	float rp, pi;
	
	dist_cylinder(p->x, q->x, rp, pi);

	hist->add(rp, pi, p->w * q->w);
      }
    }
  }
}

void corr_projected_compute_direct(Catalogues* const cats_data,
				   Catalogues* const cats_rand,
				   CorrProjected* const corr)
{
  rmax2= rp_max*rp_max + pi_max*pi_max;
  allocate_vcorr(cats_data->size());

  Histogram2D<LogBin, LinearBin>
    dd(LogBin(rp_min, rp_max, nbin), LinearBin(0.0f, pi_max, nbin_pi)),
    dr(LogBin(rp_min, rp_max, nbin), LinearBin(0.0f, pi_max, nbin_pi)),
    rr(LogBin(rp_min, rp_max, nbin), LinearBin(0.0f, pi_max, nbin_pi));


  rr.clear();
  for(Catalogues::iterator rcat= cats_rand->begin();
      rcat != cats_rand->end(); ++rcat) {

    compute_wsum(*rcat);
    count_pairs_auto_direct(*rcat, false, &rr);
      
    rr.npairs += 0.5*((*rcat)->wsum*(*rcat)->wsum - (*rcat)->w2sum);
  }

  int icat=0;
  assert(cats_data->size() > 0);
  const int rand_cats_factor= cats_rand->size() / cats_data->size();
  
  for(Catalogues::iterator cat= cats_data->begin();
      cat != cats_data->end(); ++cat) {
    dd.clear();
    dr.clear();

    // DD
    count_pairs_auto_direct(*cat, true, &dd);

    compute_wsum(*cat);
    dd.npairs += 0.5*((*cat)->wsum*(*cat)->wsum - (*cat)->w2sum);

    // DR
    // Not ndata_cat * nrand_cat crosses, but nrand_cat cross pairs
    for(size_t icat_rand=0; icat_rand<rand_cats_factor; ++icat_rand) {
      Catalogues::iterator rcat= cats_rand->begin() + icat +
	                           icat_rand * cats_data->size();
      assert(icat + icat_rand * cats_data->size() < cats_rand->size());

      if(!(*cat)->empty() && !(*rcat)->empty()) {
	count_pairs_cross_direct(*cat, *rcat, &dr);
      
	dr.npairs += (*cat)->wsum * (*rcat)->wsum;
      }
    }

    accumulate_hist(&dd);
    accumulate_hist(&dr);

    compute_corr_from_histogram2d(&dd, &dr, &rr, pi_max, vcorr.at(icat));

    icat++;
  }
  corr_projected_summarise(corr);
}


void corr_projected_compute_pairs_rr_direct(Catalogues* const cats_rand,
				   Histogram2D<LogBin, LinearBin>* const rr)
{
  msg_printf(msg_verbose, "corr_projected_compute_pairs_rr_direct\n");
  msg_printf(msg_debug, "rp_min= %e, pi_max= %e\n", rp_min, rp_max);
  rmax2= rp_max*rp_max + pi_max*pi_max;

  //
  // Count pairs RR
  //
  rr->clear();
  for(Catalogues::iterator cat=
	cats_rand->begin(); cat != cats_rand->end(); ++cat) {    
    if(!(*cat)->empty()) {
      count_pairs_auto_direct(*cat, false, rr);
      
      compute_wsum(*cat);
      rr->npairs += 0.5*((*cat)->wsum*(*cat)->wsum - (*cat)->w2sum);
    }
  }

  msg_printf(msg_verbose, "%.1lf RR pairs.\n", rr->npairs);

  accumulate_hist(rr);
}

void corr_projected_compute_pairs_all_direct(Catalogues* const cats_data,
					     Catalogues* const cats_rand,
				     Histogram2D<LogBin, LinearBin>* const dd,
				     Histogram2D<LogBin, LinearBin>* const dr,
				     Histogram2D<LogBin, LinearBin>* const rr)
{
  msg_printf(msg_verbose, "Count DD DR RR pairs (direct).\n");
  msg_printf(msg_debug, "rp_min= %e, pi_max= %e\n", rp_min, rp_max);
  rmax2= rp_max*rp_max + pi_max*pi_max;

  //
  // Count pairs RR
  //
  rr->clear();
  for(Catalogues::iterator cat=
	cats_rand->begin(); cat != cats_rand->end(); ++cat) {    
    if(!(*cat)->empty()) {
      count_pairs_auto_direct(*cat, false, rr);
      
      compute_wsum(*cat);
      rr->npairs += 0.5*((*cat)->wsum*(*cat)->wsum - (*cat)->w2sum);
    }
  }

  dd->clear();
  dr->clear();

  for(Catalogues::iterator cat= cats_data->begin();
      cat != cats_data->end(); ++cat) {
    // DD
    count_pairs_auto_direct(*cat, true, dd);

    compute_wsum(*cat);
    dd->npairs += 0.5*((*cat)->wsum*(*cat)->wsum - (*cat)->w2sum);

    // DR
    // ndata_cat * nrand_cat cross pairs
    for(Catalogues::iterator rcat= cats_rand->begin();
      rcat != cats_rand->end(); ++rcat) {    

      count_pairs_cross_direct(*cat, *rcat, dr);
      dr->npairs += (*cat)->wsum * (*rcat)->wsum;
    }
  }

  accumulate_hist(dd);
  accumulate_hist(dr);
  accumulate_hist(rr);
}

