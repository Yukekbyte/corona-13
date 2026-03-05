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
  double w_sum; // sum of weights
  double c;     // confidence weight of output (= the amount of samples behind the output sample, but can be capped)
  double W;     // contribution weight: estimate for 1/p
}
reservoir_t;

typedef struct sampler_t {
  reservoir_t **reservoirs;
}
sampler_t;

// returns random double between 0 and 1
// separate from pointsampler for things that should always be random (never stratified)
static double random_uniform() {
  return ((double)rand()/(double)RAND_MAX);
}

static void get_pixel_linear(uint64_t index, uint64_t *i, uint64_t *j, float *pixel_i, float *pixel_j) {
  pointsampler_pixel_linear(index, i, j, pixel_i, pixel_j);
}

// Updates reservoir with sample and weight.
static void update(reservoir_t *r, path_t *path, double weight, double c) {
  r->w_sum += weight;
  r->c += c;
  if (random_uniform() < weight / r->w_sum)
    path_copy(r->path, path);
}

// path must be initialized!
static md_t f(path_t *path) {
  if(is_null(path))
    return md_set1(0.0);
  return path_measurement_contribution_dx(path, 0, path->length-1);
}

// Use the integrand f as target function p_hat
// path must be initialized!
static double p_hat(path_t *path) {
  if(is_null(path))
    return 0.0;
  return md_hsum(path_measurement_contribution_dx(path, 0, path->length-1));
}

sampler_t *sampler_init() {
  uint64_t i, j;
  uint64_t w = view_width();
  uint64_t h = view_height();

  sampler_t *s = (sampler_t *)malloc(sizeof(sampler_t));
  if(s == NULL) goto fail;

  s->reservoirs = (reservoir_t**)malloc(w * sizeof(reservoir_t*));
  if(s->reservoirs == NULL) goto fail;

  for(i = 0; i < w; i++) {
    if(posix_memalign((void**)&s->reservoirs[i], 32, h * sizeof(reservoir_t)))
      goto fail;
  }

  reservoir_t *r;
  for(i = 0; i < w; i++) {
    for(j = 0; j < h; j++) {
      r = &s->reservoirs[i][j];
      if(posix_memalign((void**)&r->path, 32, sizeof(path_t))) 
        goto fail;
      set_null(r->path);
      r->w_sum = 0.;
      r->c = 0.;
      r->W = 0.;
    }
  }

  return s;

  fail:
    fprintf(stderr, "Memory allocation for reservoirs failed\n");
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
void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

static inline mf_t path_pdf_hero(const path_t *p)
{
  // this is just the hero wavelength weight:
  md_t pdf = md_set1(1.0);
  for(int v=1;v<p->length;v++)
    pdf = md_mul(pdf, mf_2d(p->v[v].pdf));

  return mf_div(md_2f(pdf), mf_set1(mf_hsum(md_2f(pdf))));
}

// Perform Resampled Importance Sampling (streaming RIS)
// r is an unused reservoir (will be reset).
static void ris(reservoir_t *r, uint64_t index, float _i, float _j) {
  // reset
  set_null(r->path);
  r->c = 0.;
  r->w_sum = 0.;
  r->W = 0.;

  const int M = 8;
  for(int i = 0; i < M; i++) {
    path_t path;
    path_init(&path, index, 0);
    path_set_pixel(&path, _i, _j);

    if(path_extend(&path)) return;

    while(1) {
      // sample light source
      if(nee_sample(&path)) break; // breaks when envmap is hit or path becomes too long
      
      // path.throughput == p_hat/pdf and mis weight = 1
      double phat = p_hat(&path);
      double pdf = md(path_pdf(&path), 0);
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
  
  // // inital candidate generation
  ris(r, path->index, _i, _j);

  // don't splat null sample
  if(is_null(r->path)) return;

  // estimator f(r.Y) * r.W
  const md_t estimator = md_mul(f(r->path), md_set1(r->W));
  const mf_t w = path_pdf_hero(r->path);

  pointsampler_splat(r->path, mf_mul(w, md_2f(estimator)));
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : ris\n");
}
