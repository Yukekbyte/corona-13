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

#define set_null(path) ((path)->length = -1)
#define is_null(path) ((path)->length == -1)
#define not_null(path) ((path)->length != -1)

// resampled importance sampling

typedef struct reservoir_t {
  path_t *path; // output sample
  md_t w_sum;   // sum of weights
  double c;     // confidence weight of output (= the amount of samples behind the output sample, but can be capped)
  md_t W;       // contribution weight: estimate for 1/p
} reservoir_t;

typedef struct sampler_t {
  reservoir_t **reservoirs;
} sampler_t;

// returns random double between 0 and 1
// separate from pointsampler for things that should always be random (never stratified)
static double random_uniform() {
  return ((double)rand()/(double)RAND_MAX);
}

static void get_pixel_linear(uint64_t index, uint64_t *i, uint64_t *j, float *pixel_i, float *pixel_j) {
  pointsampler_pixel_linear(index, i, j, pixel_i, pixel_j);
}

// Updates reservoir with sample and weight.
// Returns 1 if sample was replaced, 0 if not.
static int update(reservoir_t *r, path_t *path, md_t weight, double c) {
  r->w_sum += weight;
  r->c += c;
  if (random_uniform() < weight / r->w_sum) {
    path_copy(r->path, path);
    return 1;
  }
  return 0;
}

// path must be initialized!
static md_t f(path_t *path) {
  if(is_null(path))
    return 0.0;
  return path_measurement_contribution_dx(path, 0, path->length-1);
}

// Use the integrand f as target function p_hat
// path must be initialized!
static md_t p_hat(path_t *path) {
  if(is_null(path))
    return 0.0;
  return path_measurement_contribution_dx(path, 0, path->length-1);
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
void sampler_prepare_sample(uint64_t index) {}
void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

// Perform Resampled Importance Sampling (streaming RIS)
// r is an unused reservoir (will be reset).
static void ris(reservoir_t *r, uint64_t index, float _i, float _j) {
  // reset
  set_null(r->path);
  r->c = 0.;
  r->w_sum = 0.;
  r->W = 0.;

  const int max_length = 10;
  const int M = 8;
  for(int i = 0; i < M; i++) {
    path_t path;
    path_init(&path, index, 0);
    path_set_pixel(&path, _i, _j);

    if(path_extend(&path)) return;

    while(path.length < max_length) {
      // sample light source
      if(nee_sample(&path)) break;
      
      // path.throughput == p_hat/pdf and mis weight = 1
      md_t phat = p_hat(&path);
      md_t pdf = path_pdf(&path);
      if(phat > 0. && pdf > 0.)
        update(r, &path, phat/pdf, 1.);
      else
        r->c += 1.;
      path_pop(&path);
      
      // extend path
      if(path_extend(&path)) break;  
    }
  }

  r->w_sum *= 1./M;

  // update contribution weight W (= estimator for 1/p(r.Y)), only fails if all M samples were 0 samples
  if(not_null(r->path)) {
    r->W = r->w_sum / p_hat(r->path);
  }
}

void sampler_create_path(path_t *path)
{  
  // get pixel & reservoir from path index
  uint64_t i, j;
  float _i, _j;
  reservoir_t *r;
  get_pixel_linear(path->index, &i, &j, &_i, &_j);
  r = &rt.sampler->reservoirs[i][j];
  
  // inital candidate generation
  ris(r, path->index, _i, _j);

  // don't splat null sample
  if(is_null(r->path)) return;

  // estimator f(r.Y) * r.W
  pointsampler_splat(r->path, md_2f(f(r->path) * r->W));
}

mf_t sampler_throughput(path_t *path)
{
  // if path length is 0, return 0
  if(path->length < 1) return mf_set1(0.0f);

  // measurement contribution f (in vertex area measure)
  const md_t measurement = path_measurement_contribution_dx(path, 0, path->length-1);
  if(mf_all(mf_lte(md_2f(measurement), mf_set1(0.0f))))
    return mf_set1(0.0f);
  
  // accumulate pdf over path
  md_t pdf = md_set1(1.0);
  for(int k=0;k<path->length;k++)
    pdf = md_mul(pdf, mf_2d(path_pdf_extend(path, k)));

  // return estimate f/pdf
  return md_2f(md_div(measurement, pdf));
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
  fprintf(fd, "sampler  : ris\n");
}
