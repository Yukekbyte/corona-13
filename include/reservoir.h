#pragma once

#include "pathspace.h"
#include "points.h"

#define M 32

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
  //static const double w_arr[8] = { 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2 };
  //md_t w = md_loadu(w_arr);
  //return md_hsum(md_mul(f, w));
  return md_hsum(f);
}

float shift(path_t *shifted, pixel_t q, const path_t *source_path) {
  if(is_null(source_path)) {
    set_null(shifted);
    return 0.0;
  }

  float J = path_shift(shifted, q._i, q._j, source_path, 1);
  
  // check if shift failed
  if (J == 0.0f || p_hat(shifted) == 0.0f) {
    set_null(shifted);
    return 0.0;
  }

  return J;
}

double p_hat_from(path_t *y, pixel_t q) {
  path_t x;
  float J = shift(&x, q, y);
  if(not_null(&x))
    return p_hat(&x) * J;
  return 0;
}

double p_hat_from_opt(path_t *x, float J) {
  return p_hat(x) / J;
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
    path_init(&path, 0, 0);
    path_set_pixel(&path, q._i, q._j);
    
    if(path_extend(&path)) break;

    if(path.v[path.length-1].flags & s_environment) {
      // indicator
      r->envmap = 1;
    }
    
    // generate path tree
    while(1) {
      // sample light source
      if(nee_sample(&path)) break; // breaks when envmap is hit or path becomes too long

      // Cached value is flawed... differences give black spot artifacts
      // // use cached value instead of calculating then seperately
      // double hero_throughput = mf_hsum(path.throughput);
      // if(hero_throughput > 0.) {
      //   // w = mis * phat * 1/pdf
      //   // - mis weight is 1 for samples of same path tree (cuz each sample diff length)
      //   // - cached throughput is f / pdf 
      //   // - we use phat := hsum(f)
      //   update(r, &path, hero_throughput, 1.);
      // }
      double phat = p_hat(&path);
      double pdf = md(path_pdf(&path), 0);
      if(phat > 0. && pdf > 0.)
        update(r, &path, phat/pdf, 1.);
      else
        r->c += 1.;
      
      // splat sample anyway
      if(splat_cb) splat_cb(&path, path.throughput);
      
      path_pop(&path);
      
      // extend path
      if(path_extend(&path)) break;  
    }
  }

  // mis weights for each chosen sample per path tree is 1/M (
  // added now instead of before for simplicity, is mathematically equivalent
  r->w_sum *= 1./M;

  // update contribution weight W (= estimator for 1/p(r.Y)), only fails if all M samples were 0 samples
  if(not_null(r->path)) {
    r->W = r->w_sum / p_hat(r->path);
  } 
}
