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

// set pixel of paths ourself to control/disable anti-aliasing (necessary because we use only one reservoir per pixel)
static void get_pixel(uint64_t index, uint64_t *i, uint64_t *j) {
  // standard loop through screen pixels:
  //           w
  //       0 1 2  3  4
  //      ------------
  //    0| 0 3 6 9  12
  //  h 1| 1 4 7 10 13
  //    2| 2 5 8 11 14 (15 wraps back to index 0)
  uint64_t w = view_width();
  uint64_t h = view_height();
  *i = (index / h) % w;
  *j = index % h;
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
void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

// Perform Resampled Importance Sampling (streaming RIS)
// r is an unused reservoir (will be reset).
// init_path is an initial path for DI that has only 2 vertices starting from camera, no light source yet.
static void ris(reservoir_t *r, const path_t *init_path) {
  assert(init_path->length == 2); // only camera vertex and hitpoint
  assert(!(init_path->v[2].flags & s_environment));

  // reset
  set_null(r->path);
  r->c = 0.;
  r->w_sum = 0.;
  r->W = 0.;

  const int M = 8;
  path_t path;
  for(int k = 0; k < M; k++) {
    path_copy(&path, init_path);
    md_t pdf = 1;

    // randomly choose between one, two or three bounces
    double double_bounce_pdf = 0.33;
    double triple_bounce_pdf = 0.67;
    double u = random_uniform();
    if(double_bounce_pdf < u) {
      if(path_extend(&path)) {
        r->c += 1.;
        continue;
      }
      if(triple_bounce_pdf < u) {
        if(path_extend(&path)) {
          r->c += 1.;
          continue;
        }
        pdf *= 1 - triple_bounce_pdf;
      } else {
        pdf *= triple_bounce_pdf - double_bounce_pdf;
      }
    } else {
      pdf *= 1 - double_bounce_pdf;
    }

    // sample light source
    if(nee_sample(&path)) {
      r->c += 1.;
      continue;
    }
    
    // fast update when p_hat (f) is zero
    if(path.v[path.length-1].throughput <= 0.0) {
      r->c += 1.;
      if(path.length > 2) path_pop(&path);
      continue;
    }
      
    md_t w = 0.0;
    //md_t f = p_hat(&path); //p_hat(&path);
    //pdf *= path_pdf(&path);
    if(pdf > 0.0)
      w = (1.0/M) * (path.throughput/pdf); // 1/M uniform weights

    update(r, &path, w, 1.); // new independent sample gets confidence = 1
  }

  // update contribution weight W (= estimator for 1/p(r.Y)), only fails if all M samples were 0 samples
  if(not_null(r->path)) {
    r->W = r->w_sum / p_hat(r->path);
  }
}

void sampler_create_path(path_t *path)
{  
  // get pixel & reservoir from path index
  uint64_t i, j;
  reservoir_t *r;
  get_pixel(path->index, &i, &j);
  r = &rt.sampler->reservoirs[i][j];
  
  // extend path once to determine pixel on camera and first vertex
  path_init(path, path->index, path->sensor.camid);
  path_set_pixel(path, (float)i+0.5f, (float)j+0.5f); // +0.5 for center of pixel (no anti-aliasing!)
  if(path_extend(path)) return;

  // check for env map hit
  if(path->v[path->length-1].flags & s_environment) {
    pointsampler_splat(path, path_throughput(path));
    return;
  }

  // inital candidate generation
  ris(r, path);

  // don't splat null sample
  if(is_null(r->path)) {
    return; 
  }

  // estimator f(r.Y) * r.W
  pointsampler_splat(r->path, md_2f(f(r->path) * r->W));
  
  path_copy(path, r->path);
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
