#pragma once

#include "reservoir.h"

#define NEIGHBOUR_COUNT 4
#define NEIGHBOUR_RADIUS 10 // radius must be sufficiently big for the neighbour count.
#define PAIRWISE_COMBINE 1
#define LAMBDA_OFFSET 123.0f // is a float
#define CONFIDENCE_CAP 256. // is a double

float shift(path_t *shifted, pixel_t q, path_t *source_path) {
  if(is_null(source_path)) {
    set_null(shifted);
    return 0.0;
  }

  float J = path_shift(shifted, q._i, q._j, source_path);
  
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
  if(J > 0.)
    return p_hat(&x) * J;
  return 0.0;
}

double p_hat_from_opt(path_t *x, float J) {
  return p_hat(x) / J;
}

double mis(path_t *x, path_t *y, float J, double cy, pixel_t q[], double c[]) {
  if(is_null(y)) return 0.0;
  
  double num = cy * p_hat_from_opt(x, J);
  if(num <= 0.0) return 0.0;
  assert(num > 0.);
  
  double denom = num;
  for(int i = 0; i < NEIGHBOUR_COUNT; i++)
    denom += c[i] * p_hat_from(y, q[i]);

  if(denom <= 0.) return 0.0; // p_hat(x) should be > 0, but occasionally p_hat can be unstable, so this prevents nan's in outlier cases.

  return num / denom;
}

#if PAIRWISE_COMBINE
  // Combine, but with only two reservoirs
  static void combine_pair(pixel_t qs, reservoir_t *s, pixel_t qr, const reservoir_t *r) {
    // don't combine reservoir with itself (happens when random_neighbor fails and returns its own reservoir)
    if(s == r) { printf("tried to combine reservoir with itself\n"); return; } 

    // reservoirs can't be empty (although the path in the reservoir can still be a null sample)
    if(s->c <= 0. && r->c <= 0.) { printf("tried to combine empty reservoirs\n"); return; }
    
    path_t r_path_from_qs;
    path_t s_path_from_qr;
    double mis_r;
    double mis_s;
    float Jr = shift(&r_path_from_qs, qs, r->path);
    float Js = shift(&s_path_from_qr, qr, s->path);
    double phat_s = p_hat(s->path);
    double phat_r = p_hat(r->path);
    double phat_r_from_s = p_hat(&r_path_from_qs);
    double phat_s_from_r = p_hat(&s_path_from_qr);

    // MIS weights
    // > normally phat_s > 0. when not null (a sample can only be held if phat > 0...),
    //   but phat can occasionally evaluate to a different value as before, so to not get nan's, we also check phat_r/s <= 0.
    if(is_null(s->path) || phat_s <= 0.)
      mis_s = 0.0f;
    else
      mis_s = s->c * phat_s / (s->c * phat_s + r->c * phat_s_from_r * Js);
    
    if(is_null(r->path) || phat_r <= 0.)
      mis_r = 0.0f;
    else
      mis_r = r->c * phat_r / (r->c * phat_r + s->c * phat_r_from_s * Jr);
    
    // resampling weights
    double w_r = mis_r * phat_r_from_s * r->W * Jr;
    double w_s = mis_s * phat_s        * s->W;

    // combine reservoirs in s
    s->w_sum = w_s;
    update(s, &r_path_from_qs, w_r, r->c);

    // update estimator
    if(not_null(s->path))
      s->W = s->w_sum / p_hat(s->path);
    else
      s->W = 0.;
    
    // confidence capping
    if(s->c > CONFIDENCE_CAP) s->c = CONFIDENCE_CAP;
  }
#else
  // Combine multiple reservoirs in s
  // The samples in r[] (r[i]->path) will be shifted to the domain (pixel) of the sample in s (s->path).
  static void combine(pixel_t qs, reservoir_t *s, pixel_t qr[], reservoir_t *r[]) {
    for(int i = 0; i < NEIGHBOUR_COUNT; i++) {
      // don't combine reservoir with itself (happens when random_neighbor fails and returns its own reservoir)
      if(s == r[i]) { printf("tried to combine reservoir with itself\n"); return; }
      // reservoirs can't be empty (although the path in the reservoir can still be a null sample)
      if(s->c <= 0. && r[i]->c <= 0.) { printf("tried to combine empty reservoirs\n"); return; }
    }

    path_t Y[NEIGHBOUR_COUNT];
    float J[NEIGHBOUR_COUNT];
    double c[NEIGHBOUR_COUNT];
    double m[NEIGHBOUR_COUNT];
    double m_s;
    for(int i = 0; i < NEIGHBOUR_COUNT; i++) {
      J[i] = shift(&Y[i], qs, r[i]->path);
      c[i] = r[i]->c;
    }

    // MIS weights
    m_s = mis(s->path, s->path, 1.0f, s->c, qr, c);

    for(int i = 0; i < NEIGHBOUR_COUNT; i++) {
      pixel_t q = qr[i];
      double cy = c[i];

      qr[i] = qs;
      c[i] = s->c;
      
      m[i] = mis(r[i]->path, &Y[i], J[i], cy, qr, c);

      qr[i] = q;
      c[i] = cy;
    }

    double w[NEIGHBOUR_COUNT];
    double w_s = m_s * p_hat(s->path) * s->W;
    s->w_sum = w_s;
    for(int i = 0; i < NEIGHBOUR_COUNT; i++) {
      // resampling weight
      w[i] = m[i] * p_hat(&Y[i]) * r[i]->W * J[i];

      // combine reservoirs in s
      update(s, &Y[i], w[i], c[i]);
    }

    // update estimator
    if(not_null(s->path))
      s->W = s->w_sum / p_hat(s->path);
    else
      s->W = 0.;
    
    // confidence capping
    if(s->c > CONFIDENCE_CAP) s->c = CONFIDENCE_CAP;
  }
#endif


#if TEMPORAL_REUSE
  float shift_lambda(path_t *path) {
    if(is_null(path)) return 0.0;

    const float range = spectrum_sample_max - spectrum_sample_min;

    // shift lambdas with offset and wrap around the wavelengths that exceed spectrum_sample_max
    mf_t new_lambda = mf_add(path->lambda, mf_set1(LAMBDA_OFFSET));
    mf_t mask = mf_gt(new_lambda, mf_set1(spectrum_sample_max));
    new_lambda = mf_select(mf_sub(new_lambda, mf_set1(range)), new_lambda, mask);

    float J = path_shift_lambda(path, new_lambda);

    if(J == 0.0) set_null(path);

    return J;
  }

  // Combine reservoir s with r.
  // Assumes s and r come from the same domain (pixel)
  // Spectral shifts the sample in s.
  static void combine_temporal(reservoir_t *s, const reservoir_t *r) {
    // don't combine reservoir with itself (happens when random_neighbor fails and returns its own reservoir)
    if(s == r) { printf("tried to combine reservoir with itself (temporal)\n"); return; } 

    // reservoirs can't be empty (although the path in the reservoir can still be a null sample)
    if(s->c <= 0. && r->c <= 0.) { printf("tried to combine empty reservoirs (temporal)\n"); return; }
    
    // Shift lambda of s, should be deterministic
    float Js = shift_lambda(s->path);

    // MIS weights
    double cphats = s->c * p_hat(s->path) * Js;
    double cphatr = r->c * p_hat(r->path);
    double total = cphats + cphatr;

    if(total == 0) {
      // fast exit
      s->c += r->c;
      return;
    }

    double mis_s = cphats / total;
    double mis_r = cphatr / total;
    
    // resampling weights
    double w_r = mis_r * p_hat(r->path) * r->W;
    double w_s = mis_s * p_hat(s->path) * s->W;

    // combine reservoirs in s
    s->w_sum = w_s;
    update(s, r->path, w_r, r->c);

    // update estimator
    if(not_null(s->path))
      s->W = s->w_sum / p_hat(s->path);
    else
      s->W = 0.;
    
    // confidence capping
    if(s->c > CONFIDENCE_CAP) s->c = CONFIDENCE_CAP;
  }
#endif


