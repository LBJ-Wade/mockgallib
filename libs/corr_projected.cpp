#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "msg.h"
#include "corr_projected.h"
#include "hist2d.h"

using namespace std;

static float rmax2;
static float rp_min, rp_max, nbin, pi_max, nbin_pi;
static KDTree* tree_alloc= 0;

void corr_projected_init(const float rp_min_, const float rp_max_, const int nbin_, const float pi_max_, const int pi_nbin_)
{
  rp_min= rp_min_;
  rp_max= rp_max_;
  nbin=   nbin_;
  pi_max= pi_max_;
  nbin_pi= pi_nbin_;
}

void corr_projected_free()
{
  if(tree_alloc)
    free(tree_alloc);
}
  
//
// Local function declairations
//
static size_t count_pairs_auto(KDTree const * const tree,
			const size_t ntree,
			Histogram2D<LogBin, LinearBin>& hist);

static size_t count_pairs_cross(KDTree const * const tree1, const size_t ntree1,
			 KDTree const * const tree2,
			 Histogram2D<LogBin, LinearBin>& hist);

static void compute_corr(const Histogram2D<LogBin, LinearBin>& dd,
			 const double npairs_dd,
			 const Histogram2D<LogBin, LinearBin>& dr,
			 const double npairs_dr,
			 const Histogram2D<LogBin, LinearBin>& rr,
			 const double npairs_rr,
			 const double pi_max,
			 CorrProjected* const corr);

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
    float err=(r2 - pi*pi)/r2;
    assert(fabs(err) < 1.0e-5);
    rp= 0.0f;
  }
  else
    rp= sqrt(r2 - pi*pi);
}



size_t count_num_points(const vector<Catalogue>& v)
{
  size_t n=0;
  
  for(vector<Catalogue>::const_iterator cat= v.begin(); cat != v.end(); ++cat) {
    n += cat->size();
  }
  return n;
}

void corr_projected_compute(vector<Catalogue>& vdata,
			    vector<Catalogue>& vrandom,
			    vector<CorrProjected*>& vcorr)
{
  rmax2= rp_max*rp_max + pi_max*pi_max;

  //
  // Setup KDTree
  //
  const int nD= vdata.size(); // number of Data catalogues
  const int nR= vrandom.size(); // number of Random catalgues

  size_t nalloc= count_num_points(vdata) + count_num_points(vrandom);
  const int quota = 32;

  if(tree_alloc == 0) {
    tree_alloc= (KDTree*) malloc(sizeof(KDTree)*nalloc);
    msg_printf(msg_verbose, "%lu trees allocated (%lu Mbytes)",
	       nalloc, nalloc*sizeof(KDTree)/(1024*1024));
  }

  KDTree* tree_free= tree_alloc;
  size_t ntree_used= 0;

  // KDTree for Data
  for(vector<Catalogue>::iterator cat= vdata.begin();
      cat != vdata.end(); ++cat) {
    cat->tree= tree_free;

    cat->ntree = kdtree_construct(cat->tree, &(cat->front()),
				  cat->size(), quota);
    ntree_used += cat->ntree;
    tree_free += cat->ntree;
  }

  // KDTree for Randoms
  for(vector<Catalogue>::iterator cat= vrandom.begin();
      cat != vrandom.end(); ++cat) {
    cat->tree= tree_free;
    cat->ntree = kdtree_construct(cat->tree, &(cat->front()),
				  cat->size(), quota);
    ntree_used += cat->ntree;
    tree_free += cat->ntree;
  }

  msg_printf(msg_verbose, "%lu trees used (%lu Mbytes)",
	     ntree_used, ntree_used*sizeof(KDTree)/(1024*1024));

  
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

  double npairs_DD= 0;
  double npairs_RR= 0;
  double npairs_DR= 0;
  
  // RR
  rr.clear();
  for(vector<Catalogue>::iterator rcat= vrandom.begin();
      rcat != vrandom.end(); ++rcat) {    
    if(!rcat->empty())
      count_pairs_auto(rcat->tree, rcat->ntree, rr);

    npairs_RR += 0.5*rcat->size()*(rcat->size()-1);
  }

  int icat=0;
  for(vector<Catalogue>::iterator cat= vdata.begin();
      cat != vdata.end(); ++cat) {
    npairs_DD= 0.0;
    npairs_DR= 0.0;
    dd.clear();
    dr.clear();

    // DD
    if(!cat->empty())
      count_pairs_auto(cat->tree, cat->ntree, dd);

    npairs_DD += 0.5*cat->size()*(cat->size()-1);

    // DR
    for(vector<Catalogue>::iterator rcat= vrandom.begin();
      rcat != vrandom.end(); ++rcat) {    

      if(!cat->empty() && !rcat->empty())
	count_pairs_cross(cat->tree, cat->ntree, rcat->tree, dr);
      npairs_DR += cat->size()*rcat->size();
    }

    compute_corr(dd, npairs_DD,
		 dr, npairs_DR,
		 rr, npairs_RR,		 
		 pi_max, vcorr[icat]);

    icat++;
  }
}



static size_t count_pairs_leaf_tree_auto(KDTree const * const leaf,
					 KDTree const * const tree,
					 Histogram2D<LogBin, LinearBin>& hist)
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

    float rp, pi;

    if(leaf == tree) {
      for(Particle const *p= leaf->particles[0]; p != leaf->particles[1]; ++p) {
	for(Particle const *q= p+1; q != leaf->particles[1]; ++q) {
	  dist_cylinder(p->x, q->x, rp, pi);
	  
	  count++;
	  hist.add(rp, pi, p->w * p->w);
	}
      }
    }
    else {
      for(Particle const *p= leaf->particles[0]; p != leaf->particles[1]; ++p){
	for(Particle const *q= tree->particles[0]; q != tree->particles[1];++q){
	  dist_cylinder(p->x, q->x, rp, pi);
	  
	  count++;
	  hist.add(rp, pi, p->w * q->w);
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
			Histogram2D<LogBin, LinearBin>& hist)
{
  // Run count_pairs_leaf_tree for each leaf
  KDTree const * leaf= tree;

  size_t count= 0;

  for(size_t i=0; i<ntree; ++i) {
    if(leaf->subtree[0] == 0 && leaf->subtree[1] == 0) {
      count += count_pairs_leaf_tree_auto(leaf, tree, hist);
    }
    
    leaf++;
  }

  return count;
}

static size_t count_pairs_leaf_tree_cross(KDTree const * const leaf,
				   KDTree const * const tree,
				   Histogram2D<LogBin, LinearBin>& hist)
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

    float rp, pi;

    for(Particle const *p= leaf->particles[0]; p != leaf->particles[1]; ++p){
      for(Particle const *q= tree->particles[0]; q != tree->particles[1];++q){
	dist_cylinder(p->x, q->x, rp, pi);
	
	count++;
	hist.add(rp, pi, p->w * q->w);
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
  
size_t count_pairs_cross(KDTree const * const tree1, const size_t ntree1,
			 KDTree const * const tree2,
			 Histogram2D<LogBin, LinearBin>& hist)
{
  // Run count_pairs_leaf_tree for each leaf
  KDTree const * leaf= tree1;

  size_t count= 0;

  for(size_t i=0; i<ntree1; ++i) {
    if(leaf->subtree[0] == 0 && leaf->subtree[1] == 0) {
      count += count_pairs_leaf_tree_cross(leaf, tree2, hist);
    }
    
    leaf++;
  }

  return count;
}


// Corr Projected
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

void corr_alloc(const int n_data_cat, const int nbin, vector<CorrProjected*>& vcorr)
{
  assert(vcorr.empty());

  for(int i=0; i<n_data_cat; ++i) {
    CorrProjected* corr= new CorrProjected(nbin);
    vcorr.push_back(corr);
  }
}

void compute_corr(const Histogram2D<LogBin, LinearBin>& dd,
		  const double npairs_dd,
		  const Histogram2D<LogBin, LinearBin>& dr,
		  const double npairs_dr,
		  const Histogram2D<LogBin, LinearBin>& rr,
		  const double npairs_rr,
		  const double pi_max,
		  CorrProjected* const corr)

{
  // Compute xi(rp, pi) and project to wp(rp)
  const int nx= dd.x_nbin();
  const int ny= dd.y_nbin();
  const double dpi= 2.0*pi_max/ny;
 
  assert(corr->n == nx);

  for(int ix=0; ix<nx; ++ix) {
    corr->rp[ix]= dd.x_bin(ix);
    double wp= 0.0;
    for(int iy=0; iy<ny; ++iy) {
      double rrr= rr(ix, iy)/npairs_rr;
      if(rrr > 0.0)
	wp += (dd(ix, iy)/npairs_dd/rrr - 1.0)*dpi;

      // xi = (DD - 2*DR + RR)/RR
      // wp = \int_-pi-max^pi-max xi(rp, pi) dpi
    }    
    corr->wp[ix]= wp;
  }
}

void corr_projected_write(const int index, const vector<CorrProjected*>& vcorr)
{
  assert(!vcorr.empty());
  const int nbin= vcorr.front()->n;
  const int ndat= vcorr.size();

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

    double wp= wp_sum/ndat;
    double dwp= ndat > 1 ? sqrt((wp2_sum - ndat*wp*wp)/(ndat-1)) : 0.0;

    fprintf(fp, "%e %e %e\n", rp, wp, dwp);
  }
  fclose(fp);
}

void corr_projected_summarise(const vector<CorrProjected*>& vcorr,
			      CorrProjected* corr_out)
{
  assert(!vcorr.empty());
  const int nbin= vcorr.front()->n;
  const int ndat= vcorr.size();

  for(int i=0; i<nbin; ++i) {
    double rp= vcorr.front()->rp[i];
    double wp_sum= 0.0;
    double wp2_sum= 0.0;
    for(vector<CorrProjected*>::const_iterator corr= vcorr.begin();
	corr != vcorr.end(); ++corr) {
      wp_sum +=  (*corr)->wp[i];
      wp2_sum += (*corr)->wp[i]*(*corr)->wp[i];
    }

    double wp= wp_sum/ndat;
    double dwp= ndat > 1 ? sqrt((wp2_sum - ndat*wp*wp)/(ndat-1)) : wp;

    corr_out->rp[i]= rp;
    corr_out->wp[i]= wp;
    corr_out->dwp[i]= dwp;
  }
}

