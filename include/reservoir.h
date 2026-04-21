#pragma once

#include "pathspace.h"
#include "points.h"

#define M 8

#define set_null(path) ((path)->length = -1)
#define is_null(path) ((path)->length == -1)
#define not_null(path) ((path)->length != -1)

typedef struct reservoir_t {
  path_t *path; // output sample
  double w_sum; // sum of weights
  double c;     // confidence weight of output (= the amount of samples behind the output sample)
  double W;     // contribution weight: estimate for 1/p
  uint8_t envmap;
} reservoir_t;

typedef struct pixel_t {
  uint64_t i;
  uint64_t j;
  float _i;
  float _j;
} pixel_t;

// Updates reservoir with sample and weight.
static void update(reservoir_t *r, path_t *path, double weight, double c) {
  r->w_sum += weight;
  r->c += c;
  if (points_rand(rt.points, common_get_threadid()) * r->w_sum < weight)
    path_copy(r->path, path);
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

typedef void (*splat_fn)(path_t*, mf_t);

// Perform Resampled Importance Sampling (streaming RIS)
// r will be reset and filled with initial samples for pixel q
static void ris(pixel_t q, reservoir_t *r, splat_fn splat_cb) {
  // reset
  set_null(r->path);
  r->c = 0.;
  r->w_sum = 0.;
  r->W = 0.;
  r->envmap = 0;

  for(int i = 0; i < M; i++) {
    path_t path;

    // random path index (32 bit suffices (instead of 64)), needed for different lambda sampling
    uint32_t index = (uint32_t)(points_rand(rt.points, common_get_threadid()) * 4294967296.0f); // 2^32
    
    path_init(&path, index, 0);
    path_set_pixel(&path, q._i, q._j);
    path_set_aperture(&path, 0, 0);
    
    if(path_extend(&path)) break;

    if(path.v[path.length-1].flags & s_environment) 
      r->envmap = 1; // indicator
    
    // generate path tree
    while(1) {
      // sample light source
      if(nee_sample(&path)) break; // breaks when envmap is hit or path becomes too long

      double phat = p_hat(&path);
      double pdf = md_hsum(path_pdf(&path));
      if(phat > 0. && pdf > 0.)
        update(r, &path, phat/pdf, 1.);
      else
        r->c += 1.;
      
      // splat sample anyway
      mf_t mis = sampler_mis(&path); // Hero MIS, don't confuse this with the resampling weight MIS from above (which is 1/M)!
      if(splat_cb) splat_cb(&path, mf_mul(path.throughput, mis));
      
      path_pop(&path);
      
      // extend path
      if(path_extend(&path)) break;  
    }
  }

  // mis weights for each chosen sample per path tree is 1/M
  // added now instead of before for simplicity, because it's mathematically equivalent
  r->w_sum *= 1./M;

  // update contribution weight W (= estimator for 1/p(r.Y)), only fails if all M samples were 0 samples
  if(not_null(r->path)) {
    r->W = r->w_sum / p_hat(r->path);
  } 
}
