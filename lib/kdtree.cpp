#include <iostream>
#include <cmath>
#include <cassert>
#include "particle.h"
#include "kdtree.h"

// KDTree constructions
using namespace std;

static KDTree* root; // debug!!

static void compute_bounding_box(Particle const * const p, const index_t n, float* const left, float* const right);

static index_t split_particles(Particle* const particle,
			       const index_t particle_begin, 
			       const index_t particle_end, 
			       const float mid, const int i);

static KDTree* construct_tree_recursive(KDTree* const tree,
					Particle* const particle,
					const index_t particle_begin,
					const index_t particle_end,
					const float left[],
					const float right[],
					const int quota);

size_t kdtree_construct(KDTree* const tree, Particle* const particle, const index_t np, const int quota)
{
  if(np <= 0) return 0;
  root=tree;//debug;
  float left[3], right[3];
  compute_bounding_box(particle, np, left, right);
  
  KDTree* const tree_end=
    construct_tree_recursive(tree, particle, 0, np, left, right, quota);

  return tree_end - tree;
}

KDTree* construct_tree_recursive(KDTree* const tree,
				 Particle* const particle,
				 const index_t particle_begin,
				 const index_t particle_end,
				 const float left[],
				 const float right[],
				 const int quota)
{
  KDTree* tree_free= tree + 1;

  for(int k=0; k<3; ++k) {
    //printf("%e %e\n", left[0], right[0]);
    tree->left[k]= left[k];
    tree->right[k]= right[k];
  }

  assert(tree->right[0] - tree->left[0] > 0.0f);
  //printf("tree %lu %d\n", tree-root, particle_end - particle_begin);

  /*
  if(particle_begin - particle_end < 40) {
    for(int i=particle_begin; i<particle_end; ++i)
      printf("%e %e %e\n",
	     particle[i].x[0], particle[i].x[1], particle[i].x[2]);
  }
  */

  
  int i= 0;
  float dx= right[0] - left[0];
  for(int k=1; k<3; ++k) {
    float dx_k= right[k] - left[k];
    if(dx_k > dx) {
      i= k; dx= dx_k;
    }
  }
  //cerr << "i= " << i << endl;
  //fprintf(stderr, "%e %e %e\n", right[0] - left[0],
  //right[1] - left[1], right[2] - left[2]);

  tree->subtree[0]= tree->subtree[1]= 0;
  tree->particles[0]= particle + particle_begin;
  tree->particles[1]= particle + particle_end;
  
  // terminal condition
  if(particle_end - particle_begin < quota) {
    return tree_free;
  }

  //const float mid= 0.5f*(left[i]+right[i]);
  const float mid= left[i] + 0.5f*(right[i] - left[i]);
  assert(left[i] <= mid && mid <= right[i]);
  const index_t particle_mid=
    split_particles(particle, particle_begin, particle_end, mid, i);

  // Note: you can improve this by reusing this tree when one of the
  // subtree is empty
  
  if(particle_mid != particle_begin) {
    float right_next[]= {right[0], right[1], right[2]};
    right_next[i]= mid;
    tree->subtree[0]= tree_free++; 
    tree_free= construct_tree_recursive(tree->subtree[0], particle,
					particle_begin, particle_mid,
					left, right_next, quota);
  }
  if(particle_mid < particle_end) {
    float left_next[]= {left[0], left[1], left[2]};
    left_next[i]= mid;
    tree->subtree[1]= tree_free++; 
    tree_free= construct_tree_recursive(tree->subtree[1], particle,
					particle_mid, particle_end,
					left_next, right, quota);
  }

  // Caution; boundary not checked for tree_free++ 
  // unpredictable result (possible segmentation fault) if not enough trees
  // are allocated
  assert(tree->subtree[0] || tree->subtree[1]);

  return tree_free;
}

index_t split_particles(Particle* const particle,
			const index_t particle_begin, 
			const index_t particle_end, 
			const float mid, const int i)
{
  // Rearranges the particles into two groups
  // particle[particle_begin ... j1] with x[i] < mid
  // particle[j1 ... particle_end] with x[i] >= mid
  
  assert(particle_begin < particle_end);
  int j1= particle_begin;
  int j2= particle_end-1;
  
  while(1) {
    while(j1 < particle_end && particle[j1].x[i] < mid)
      ++j1;
    while(j2 >= particle_begin && particle[j2].x[i] >= mid)
      --j2;
    if(j1 < j2) {
      const Particle temp= particle[j1];
      particle[j1]= particle[j2];
      particle[j2]= temp;
      ++j1; --j2;
    }
    else
      break;
  }
  assert(j1 >= particle_begin && j1 <= particle_end && (j1 - j2 == 1));
  return j1;
}


void compute_bounding_box(Particle const * const p, const index_t np, float* const left, float* const right)
{
  for(int k=0; k<3; ++k) {
    left[k]= right[k]= p->x[k];
  }
  
  for(int i=0; i<np; ++i) {
    for(int k=0; k<3; ++k) {
      if(p[i].x[k] < left[k]) left[k]= p[i].x[k];
      if(p[i].x[k] > right[k]) right[k]= p[i].x[k];
    }
  }
}
