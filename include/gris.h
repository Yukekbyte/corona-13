#pragma once

#include "reservoir.h"


float shift(path_t *shifted, pixel_t q, const path_t *source_path) {
  if(is_null(source_path)) {
    set_null(shifted);
    return 0.0;
  }

  float i = q.i + fractf(source_path->sensor.pixel_i);
  float j = q.j + fractf(source_path->sensor.pixel_j);
  
  float J = path_shift(shifted, i, j, source_path);
  
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

double mis(path_t *x, path_t *y, float J, int i, pixel_t q[], double c[]) {
  if(is_null(y)) return 0.0;
  
  double num = c[i] * p_hat_from_opt(x, J);
  if(num <= 0.0) return 0.0;
  
  double denom = num;
  for(int j = 0; j < K; j++) {
    if(j == i) continue;
    denom += c[j] * p_hat_from(y, q[j]);
  }

  if(denom <= 0.) return 0.0; // paranoid case because p_hat(x) should be > 0

  return num / denom;
}

#if PAIRWISE_COMBINE
  // Combine, but with only two reservoirs
  static void combine_pair(reservoir_t *s, pixel_t qs, reservoir_t *r, pixel_t qr) {
    assert(s != r);

    // reservoirs can't be empty (although the path in the reservoir can still be a null sample)
    if(s->c <= 0. && r->c <= 0.) { printf("tried to combine empty reservoirs\n"); return; }
    
    path_t r_path_in_s, s_path_in_r;
    double mis_r, mis_s;
    float Jr = shift(&r_path_in_s, qs, r->path);
    float Js = shift(&s_path_in_r, qr, s->path);
    double phat_s = p_hat(s->path);
    double phat_r = p_hat(r->path);
    double phat_r_in_s = p_hat(&r_path_in_s);
    double phat_s_in_r = p_hat(&s_path_in_r);

    // MIS weights
    if(phat_s == 0.0)
      mis_s = 0.0;
    else
      mis_s = s->c * phat_s / (s->c * phat_s + r->c * phat_s_in_r * Js);
    
    if(phat_r == 0.0 || Jr == 0.0)
      mis_r = 0.0;
    else
      mis_r = (r->c * phat_r / Jr) / (r->c * phat_r / Jr + s->c * phat_r_in_s);
    
    // resampling weights
    double w_r = mis_r * phat_r_in_s * r->W * Jr;
    double w_s = mis_s * phat_s      * s->W;

    // combine reservoirs in s
    s->w_sum = w_s;
    update(s, &r_path_in_s, w_r, r->c);

    // update estimator
    if(not_null(s->path))
      s->W = s->w_sum / p_hat(s->path);

    // confidence capping
    if(s->c > CONFIDENCE_CAP) s->c = CONFIDENCE_CAP;
  }
#else
  // Combine multiple reservoirs in r
  // The samples in x[] (x[i]->path) will be shifted to q.
  static void combine(reservoir_t *r, pixel_t q, reservoir_t *x[], pixel_t qs[]) {
    // I suppose we could assert that the reservoirs are different here as well...

    for(int i = 0; i < K; i++) {
      // reservoirs can't be empty (although the path in the reservoir can still be a null sample)
      if(x[i]->c <= 0.) { printf("tried to combine empty reservoirs\n"); return; }
    }

    path_t Y[K];
    float J[K];
    double c[K];
    double m[K];

    for(int i = 0; i < K; i++) {
      J[i] = shift(&Y[i], q, x[i]->path);
      c[i] = x[i]->c;
    }

    // MIS weights
    for(int i = 0; i < K; i++)
      m[i] = mis(x[i]->path, &Y[i], J[i], i, qs, c);

    reset(r);
    
    // Resampling weights
    double w[K];
    for(int i = 0; i < K; i++) {
      if(not_null(&Y[i]))
        w[i] = m[i] * p_hat(&Y[i]) * x[i]->W * J[i];
      else
        w[i] = 0.0;
      update(r, &Y[i], w[i], c[i]);
    }

    // update estimator
    if(not_null(r->path))
      r->W = r->w_sum / p_hat(r->path);

    // confidence capping
    if(r->c > CONFIDENCE_CAP) r->c = CONFIDENCE_CAP;
  }
#endif


#if TEMPORAL_REUSE
  float spectral_shift(path_t *shifted, path_t *source, uint8_t invert) {
    if(is_null(source)) {
      set_null(shifted);
      return 0.0;
    }

    const float range = spectrum_sample_max - spectrum_sample_min;
    float offset = invert ? -LAMBDA_OFFSET : LAMBDA_OFFSET; //(rt.points, common_get_threadid()) * range;

    // shift lambdas with offset and wrap around the wavelengths that exceed spectrum_sample_max
    mf_t new_lambda = mf_add(source->lambda, mf_set1(offset));
    mf_t mask = mf_gt(new_lambda, mf_set1(spectrum_sample_max));
    new_lambda = mf_select(mf_sub(new_lambda, mf_set1(range)), new_lambda, mask);

    float J = path_shift_lambda(shifted, new_lambda, source);

    if(J == 0.0) set_null(shifted);

    return J;
  }

  void spectral_combine(reservoir_t *s, reservoir_t *r) {
    assert(s != r);

    if(s->c <= 0. && r->c <= 0.) { printf("Tried to combine empty reservoirs (temporal)\n"); return; }

    path_t rlambda;
    float Jr = spectral_shift(&rlambda, r->path, 0);
    double mis_r, mis_s;

    mis_r = 0.5;
    mis_s = 0.5;

    double w_r = mis_r * p_hat(&rlambda) * r->W * Jr;
    double w_s = mis_s * p_hat(s->path) * s->W;

    s->w_sum = w_s;
    update(s, &rlambda, w_r, r->c);

    if(not_null(s->path))
      s->W = s->w_sum / p_hat(s->path);
    else
      s->W = 0.;

    if(s->c > CONFIDENCE_CAP) s->c = CONFIDENCE_CAP;
  }
#endif


