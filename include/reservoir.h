#pragma once

#include "pathspace.h"
#include "pointsampler.h"
#include "points.h"

#define M 16                     // Amount of NEE-trees an intial reservoir is filled with.
#define SPECTRAL_REUSE 1        // Turns Spectral Re-use on/off
#define SPATIAL_REUSE_PASSES 4  // Amount of spatial re-use passes.

#define NEIGHBOUR_COUNT 4   // Amount of neighbours to sample and combine with during a spatial re-use pass.
#define NEIGHBOUR_RADIUS 10 // Half the length of the side of the square where neighbouring pixels are sampled in. Must be sufficiently big for the neighbour count!
#define PAIRWISE_COMBINE 1  // Combine neighbouring reservoirs in pairs instead of all at once.
#define PAIRWISE_MIS 1      // Use linear pairwise MIS weights (instead of balance heuristic = quadratic) only has effect if PAIRWISE_COMBINE is 0.
#define LAMBDA_OFFSET 20.0f // The wavelength perturbation used in temporal re-use
#define CONFIDENCE_CAP 200.  // The confidence weight limit.

#define K (NEIGHBOUR_COUNT + 1)
#define set_null(path) ((path)->length = -1)
#define is_null(path) ((path)->length == -1)
#define not_null(path) ((path)->length != -1)

typedef struct reservoir_t {
  path_t *path; // output sample
  double w_sum; // sum of weights
  double c;     // confidence weight of output (= the amount of samples behind the output sample)
  double W;     // contribution weight: estimate for 1/p
} reservoir_t;

typedef struct pixel_t {
  int i;
  int j;
} pixel_t;

// fractional part of float
static inline float fractf(float x) { return x - (int)x; }

// Updates reservoir with sample and weight.
static void update(reservoir_t *r, path_t *path, double weight, double c) {
  r->w_sum += weight;
  r->c += c;
  if (points_rand(rt.points, common_get_threadid()) * r->w_sum < weight)
    path_copy(r->path, path);
}

static inline void reset(reservoir_t *r) {
  set_null(r->path);
  r->c = 0.;
  r->w_sum = 0.;
  r->W = 0.;
}

// Use the integrand f as target function p_hat
static double p_hat(path_t *path) {
  if(is_null(path))
    return 0.0;
  md_t f = path_measurement_contribution_dx(path, 0, path->length-1);
  return md_hsum(f);
}

static inline mf_t sampler_mis(const path_t *p) 
{
  md_t pdf = md_set1(1.0);
  for(int v=1;v<p->length;v++)
    pdf = md_mul(pdf, mf_2d(p->v[v].pdf));

  return mf_div(md_2f(pdf), mf_set1(mf_hsum(md_2f(pdf))));
}

static void ris_update(reservoir_t *r, path_t *path) {
  double phat = p_hat(path);
  double pdf = md_hsum(path_pdf(path));
  if(phat > 0. && pdf > 0.)
    update(r, path, phat/pdf, 1.);
  else
    r->c += 1.;
}

// Perform Resampled Importance Sampling (streaming RIS)
// r will be reset and filled with initial samples for pixel q
static void ris(pixel_t q, reservoir_t *r) {
  reset(r);

  for(int i = 0; i < M; i++) {
    path_t path;

    // random path index (32 bit suffices (instead of 64)), needed for different lambda sampling
    uint32_t index = (uint32_t)(points_rand(rt.points, common_get_threadid()) * 4294967296.0f); // 2^32
    
    path_init(&path, index, 0);
    // anti aliasing
    float pixel_i, pixel_j;
    pointsampler_subpixel(q.i, q.j, &pixel_i, &pixel_j);
    path_set_pixel(&path, pixel_i, pixel_j);
    
    if(path_extend(&path)) assert(0);

    // for light sources as hit point
    ris_update(r, &path);
    
    // generate path tree
    while(1) {
      // sample light source
      if(nee_sample(&path)) break; // breaks when envmap is hit or path becomes too long

      ris_update(r, &path);
      
      path_pop(&path);
      
      // extend path
      if(path_extend(&path)) {
        r->c += 1.;
        break;
      }  
    }
  }

  // mis weights for each chosen sample per path tree is 1/M
  // added now instead of before for simplicity, because it's mathematically equivalent
  r->w_sum *= 1./M;

  // update contribution weight W (= estimator for 1/p(r.Y)), only fails if all M samples were 0 samples
  if(not_null(r->path)) {
    r->W = r->w_sum / p_hat(r->path);
  }

  if(r->c >= CONFIDENCE_CAP) r->c = CONFIDENCE_CAP;
}
