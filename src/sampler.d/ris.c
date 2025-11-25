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

// resampled importance sampling

typedef struct reservoir_t {
  path_t *path; // output sample
  double w_sum; // sum of weights
  double c; // confidence weight of output
  mf_t W; // contribution weight: estimate for 1/f (because we choose p_hat := f)
} reservoir_t;

typedef struct sampler_t {
  reservoir_t **reservoirs;
} sampler_t;

// Updates reservoir with sample and weight.
// Returns 1 if sample was replaced, 0 if not.
static int update(reservoir_t *r, path_t *path, double weight, double c) {
  r->w_sum += weight;
  r->c += c;
  if (((double)rand()/(double)RAND_MAX) < weight / r->w_sum) {
    path_copy(r->path, path);
    return 1;
  }
  return 0;
}

static void combine(reservoir_t *s, const reservoir_t *r1, const reservoir_t *r2) {
  s->c = 0;
  s->w_sum = 0;
  s->W = 0;

  //TODO
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
      r->w_sum = 0;
      r->c = 0;
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

static inline mf_t sampler_mis(const path_t *p)
{
  // this is just the hero wavelength weight:
  md_t pdf = md_set1(1.0);
  for(int v=1;v<p->length;v++)
    pdf = md_mul(pdf, mf_2d(p->v[v].pdf));

  return mf_div(md_2f(pdf), mf_set1(mf_hsum(md_2f(pdf))));
}

// Perform Resampled Importance Sampling (streaming RIS)
// r is an unused reservoir (will be reset).
// init_path is an initial path for DI that has only 2 vertices starting from camera, no light source yet.
static void ris(reservoir_t *r, const path_t *init_path) {
  assert(init_path->length == 2); // only camera vertex and hitpoint
  assert(!(init_path->v[0].mode & s_emit));

  // reset
  r->c = 0;
  r->w_sum = 0;
  r->W = 0;

  int M = 8;
  for(int k = 0; k < M; k++) {
    path_t p; // declaring this out of loop creates issues...
    path_copy(&p, init_path);
    
    // direct illumination, fails when hitpoint of init_path on envmap
    if(nee_sample(&p)) continue;

    double w = 1.0/M * path_throughput(&p); // 1/M uniform weights
    
    update(r, &p, w, 1); // new independent sample gets confidence = 1
  }

  // if sample exists
  if(r->w_sum > 0) {
    r->W = mf_mul(mf_div(mf_set1(1.0f), path_throughput(r->path)), mf_set1(r->w_sum));
  }

  // confidence capping
  if(r->c > 50) r->c = 50;
}

void sampler_create_path(path_t *path)
{  
  // init and extend path once to determine pixel on camera and first vertex (stuff handled by pointsampler)
  path_init(path, path->index, path->sensor.camid);
  if(path_extend(path)) return;

  // check for env map hit
  if(path->v[path->length-1].flags & s_environment) {
    pointsampler_splat(path, path_throughput(path));
    return;
  }

  reservoir_t *r;
  uint64_t i = (uint64_t)path->sensor.pixel_i;
  uint64_t j = (uint64_t)path->sensor.pixel_j;
  r = &rt.sampler->reservoirs[i][j];

  // inital candidate generation
  reservoir_t rris;
  rris.path = (path_t*)malloc(sizeof(path_t));

  ris(&rris, path);

  // weighted throughput
  pointsampler_splat(rris.path, mf_mul(path_throughput(rris.path), rris.W));
  
  path_copy(path, rris.path);

  free(rris.path);
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
