/*
    This file is part of corona-13.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "pathspace/nee.h"
#include "view.h"
#include "sampler.h"
#include "spectrum.h"
#include "pointsampler.h"
#include "pathspace/multichain.h"

#define set_null(path) ((path)->length = -1)
#define is_null(path) ((path)->length == -1)
#define not_null(path) ((path)->length != -1)

// Reservoir-based Spatio-Temporal Importance Resampling (ReSTIR)

typedef struct reservoir_t {
  path_t *path; // output sample
  md_t w_sum;   // sum of weights
  double c;     // confidence weight of output (= the amount of samples behind the output sample)
  md_t W;       // contribution weight: estimate for 1/p
} reservoir_t;

typedef struct pixel_t {
  uint64_t i;
  uint64_t j;
  float _i;
  float _j;
} pixel_t;

typedef struct sampler_t {
  reservoir_t **reservoirs;
} sampler_t;

// returns random double between 0 and 1
// separate from pointsampler for things that should always be random (never stratified)
static double random_uniform() {
  return ((double)rand()/(double)RAND_MAX);
}

// set pixel of paths ourself to: 
// - avoid race conditions between threads by using the halton sequence
// - control/disable anti-aliasing (necessary because we use only one reservoir per pixel)
static void get_pixel(const uint64_t index, pixel_t *q, const int set) {
  if(set) {
    // // linear sequence
    // uint64_t w = view_width();
    // uint64_t h = view_height();
    // q->i = (index / h) % w;
    // q->j = index % h;

    // halton sequence
    pointsampler_pixel(index, &q->_i, &q->_j);
    q->i = (uint64_t)(q->_i);
    q->j = (uint64_t)(q->_j);
  }

  // disable anti-aliasing
  q->_i = (float)(q->i) + 0.5f;
  q->_j = (float)(q->j) + 0.5f;
}

static void path_from_pixel(const pixel_t q, path_t *path) {
  path_init(path, 0, 0); // TODO potentially update camid if using other cameras (see sampler_create_path)
  path_set_pixel(path, q._i, q._j);
  path_extend(path);
}

// Updates reservoir with sample and weight.
// Returns 1 if sample was replaced, 0 if not.
static int update(reservoir_t *r, path_t *path, md_t weight, double c) {
  r->w_sum += weight;
  r->c += c;
  if (random_uniform() * r->w_sum < weight) {
    path_copy(r->path, path);
    return 1;
  }
  return 0;
}

static md_t f(path_t *path) {
  if(is_null(path))
    return 0.0;
  return path_measurement_contribution_dx(path, 0, path->length-1);
}

// Use the integrand f as target function p_hat
static md_t p_hat(path_t *path) {
  if(is_null(path))
    return 0.0;
  return path_measurement_contribution_dx(path, 0, path->length-1);
}

static md_t p_hat_from(const pixel_t q, path_t *path) {
  if(is_null(path))
    return 0.0;

  path_t path_inverted;
  float J = path_shift(&path_inverted, q._i, q._j, path, 1);
  return p_hat(&path_inverted) * J;
}

// Perform Resampled Importance Sampling (streaming RIS)
// r will be reset and filled with initial samples for pixel q
static void ris(const pixel_t q, reservoir_t *r) {
  
  // reset
  set_null(r->path);
  r->c = 0.;
  r->w_sum = 0.;
  r->W = 0.;
  
  path_t path;
  path_from_pixel(q, &path);
  
  // early return when path on envmap
  if(path.v[path.length-1].flags & s_environment) {
    // give r no confidence but manually update r->path
    // this means r->path is not a null sample but r->c is zero! (tricky case...)
    path_copy(r->path, &path);
    return;
  }
  
  const int M = 8;
  for(int k = 0; k < M; k++) {

    // direct illumination
    if(nee_sample(&path)) {
      r->c += 1.; // fast update when sampling fails
      if(path.length > 2) path_pop(&path);
      continue;
    }

    // if(!path_visible(&path, 2)) {
    //   r->c += 1.;
    //   path_pop(&path);
    //   continue;
    // }

    md_t w = 0.0;
    md_t f = p_hat(&path);
    md_t p = path_pdf(&path);
    if(p > 0.0)
      w = (1.0/M) * (f/p); // 1/M uniform weights

    update(r, &path, w, 1.); // new independent sample gets confidence = 1
    path_pop(&path);
  }

  // update contribution weight W (= estimator for 1/p(r.Y)), only fails if all M samples were 0 samples
  if(not_null(r->path)) {
    r->W = r->w_sum / p_hat(r->path);
  } 
}

float path_shift_wrapper(path_t *shifted, const pixel_t q, const path_t *source_path) {
  if(is_null(source_path)) {
    set_null(shifted);
    return 0.0;
  }

  return path_shift(shifted, q._i, q._j, source_path, 1);
}

// Combine reservoir s with r.
// The sample in r (r->path) will be shifted to the domain (pixel) of the sample in s (s->path).
static void combine(const pixel_t qs, reservoir_t *s, pixel_t const qr, const reservoir_t *r) {
  // don't combine reservoir with itself (happens when random_neighbor fails and returns its own reservoir)
  if(s == r) { printf("tried to combine reservoir with itself\n"); return; } 

  // reservoirs can't be empty (although the path in the reservoir can still be a null sample)
  if(s->c <= 0. && r->c <= 0.) { printf("tried to combine empty reservoirs\n"); return; }

  if(is_null(r->path)) {
    s->c += r->c;
    return;
  }
  
  path_t r_path_from_qs;
  md_t mis_r;
  md_t mis_s;
  md_t w_r;
  md_t w_s;

  // shift r's path to s's domain
  float J = path_shift_wrapper(&r_path_from_qs, qs, r->path);

  // option 1
  if(is_null(s->path)) {
    s->c += r->c;
    return;
  }

  // option 2
  // if(is_null(s->path)) {
  //   path_copy(s->path, &r_path_from_qs);
  //   s->w_sum = p_hat(&r_path_from_qs) * r->W * J;
  //   s->c += r->c;
  //   s->W = s->w_sum / p_hat(s->path);
  //   return;
  // }

  mis_r = (r->c) / (s->c + r->c);
  mis_s = (s->c) / (s->c + r->c);
  
  w_r = mis_r * p_hat(&r_path_from_qs) * r->W * J;
  w_s = mis_s * p_hat(s->path) * s->W;

  // combine reservoirs in s
  s->w_sum = w_s;
  update(s, &r_path_from_qs, w_r, r->c);

  // update estimator
  if(not_null(s->path))
    s->W = s->w_sum / p_hat(s->path); // NaN if s->path is null (not really a problem)
  else
    s->W = 0.;
  
  // confidence capping
  if(s->c > 100.) s->c = 100.;
}

static void random_neighbor(const pixel_t q, const path_t *path, reservoir_t **n, pixel_t *q_n) {
  uint64_t w = view_width();
  uint64_t h = view_height();

  int l, m; // int instead of uint64_t so they can be negative (important for clamping)

  const int d = 30; // sample in 2*d x 2*d square with (i, j) in the center
  const int MAX_ATTEMPTS = 10; // when i == l && j == m or neighbor not geometrically similar enough
  for(int k = 0; k < MAX_ATTEMPTS; k++) {
    // TODO: make edge clamping uniform probablility
    l = q.i + (2*d*random_uniform() - d);
    m = q.j + (2*d*random_uniform() - d);

    l = CLAMP(l, 0, w-1);
    m = CLAMP(m, 0, h-1);

    q_n->i = (uint64_t)l;
    q_n->j = (uint64_t)m;

    if(q_n->i == q.i && q_n->j == q.j) continue;

    // TODO: filter geometric similarity
    // r->path->v[1] must be similar to path->v[1]

    *n = &rt.sampler->reservoirs[q_n->i][q_n->j];
    get_pixel(0, q_n, 0); // only set _i, and _j

    return;
  }

  *n = &rt.sampler->reservoirs[q.i][q.j]; // failed: return same reservoir
  // q_n unchanged
  return;
}
    
sampler_t *sampler_init() {
  uint64_t i, j;
  uint64_t w = view_width(); // 1024 (standard)
  uint64_t h = view_height(); // 576 (standard)

  sampler_t *s = (sampler_t *)malloc(sizeof(sampler_t));

  s->reservoirs = (reservoir_t**)malloc(w * sizeof(reservoir_t*));
  for(i = 0; i < w; i++) {
    s->reservoirs[i] = (reservoir_t*)malloc(h * sizeof(reservoir_t));
  }

  reservoir_t *r;
  for(i = 0; i < w; i++) {
    for(j = 0; j < h; j++) {
      r = &s->reservoirs[i][j];
      r->path = (path_t *)malloc(sizeof(path_t));
      set_null(r->path);
      r->w_sum = 0.;
      r->c = 0.;
      r->W = 0.;
    }
  }

  return s;
}

void sampler_cleanup(sampler_t *s) {
  uint64_t i, j;
  uint64_t w = view_width();
  uint64_t h = view_height();
  
  for(i = 0; i < w; i++) {
    for(j = 0; j < h; j++) {
      free(s->reservoirs[i][j].path);
    }
    free(s->reservoirs[i]);
  }
  free(s->reservoirs);
  free(s);
}

void sampler_prepare_frame(sampler_t *s) { 
  uint64_t w = view_width();
  uint64_t h = view_height();
  uint64_t i, j;
  pixel_t q;
  
  // for each pixel on screen, perform initial RIS
  for(i = 0; i < w; i++) {
    for(j = 0; j < h; j++) {
      q.i = i;
      q.j = j;
      //q._i = (float)(i) + 0.5f;
      //q._j = (float)(j) + 0.5f;
      get_pixel(0, &q, 0); // only set _i, and _j
      ris(q, &s->reservoirs[i][j]);
    }
  }
}

void sampler_clear(sampler_t *s) {}

void sampler_create_path(path_t *path)
{  
  // get pixel & reservoir from path index
  reservoir_t *r;
  pixel_t q;
  get_pixel(path->index, &q, 1);
  r = &rt.sampler->reservoirs[q.i][q.j];
  r->path->index = path->index;
  //r->path->sensor.camid = path->sensor.camid; // TODO potentially if using other cameras (see path_from_pixel)

  if(not_null(r->path)) {
    assert(r->path->sensor.pixel_i == q._i);
    assert(r->path->sensor.pixel_j == q._j);
  }

  // generate initial candidates
  //ris(q, r); //done in sampler_prepare_frame()

  // check for env map hit
  if(not_null(r->path) && r->path->v[r->path->length-1].flags & s_environment) {
    pointsampler_splat(r->path, path_throughput(r->path));
    path_copy(path, r->path);
    return;
  }

  #if(1)
  // spatial re-use
  const int neighbors = 10;
  {
  reservoir_t neighbor;
  reservoir_t *n;
  pixel_t q_n;
  neighbor.path = (path_t *)malloc(sizeof(path_t));

  for(int k = 0; k < neighbors; k++) {
    // TODO: pick neighbors with low discrepancy sequence
    random_neighbor(q, r->path, &n, &q_n);
    neighbor.w_sum = n->w_sum;
    neighbor.c = n->c; 
    neighbor.W = n->W;
    path_copy(neighbor.path, n->path); // copy to avoid race conditions/dirty reads (still not fully bullet proof as the path can have changed between w_sum copy and now)
    
    combine(q, r, q_n, &neighbor);
  }

  free(neighbor.path);
  }
  #endif

  // don't splat null sample
  if(is_null(r->path)) return;

  // estimator f(r.Y) * r.W
  pointsampler_splat(r->path, md_2f(f(r->path) * r->W));
  
  path_copy(path, r->path);
}

mf_t sampler_throughput(path_t *path)
{
  // if path length is 0, return 0
  if(path->length < 1) 
    return 0;

  // measurement contribution f (in vertex area measure)
  const md_t measurement = f(path);
  if(measurement <= 0.)
    return 0;
  
  // accumulate pdf over path
  md_t p = 1.;
  for(int k=0; k < path->length; k++)
    p = p * path_pdf_extend(path, k);

  // return estimate f/pdf
  return measurement / p;
}

md_t sampler_mis_weight(path_t *p)
{
  return md_set1(1.0);
}

md_t sampler_sum_pdf_dwp(path_t *p)
{
  // only used by hwl, hrec guided
  return md_set1(1.0);
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : restir\n");
}
