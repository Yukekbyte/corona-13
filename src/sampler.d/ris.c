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
  md_t w_sum; // sum of weights
  double c; // confidence weight of output (= the amount of samples behind the output sample)
  md_t W; // contribution weight: estimate for 1/p
} reservoir_t;

typedef struct sampler_t {
  reservoir_t **reservoirs;
} sampler_t;

// returns random double between 0 and 1
static double random_uniform() {
  return ((double)rand()/(double)RAND_MAX);
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

static md_t f(path_t *path) {
  if(path->length != 3)
    return 0.0;
  return path_measurement_contribution_dx(path, 0, path->length-1);
}

// Use the integrand f as target function p_hat
static md_t p_hat(path_t *path) {
  if(path->length != 3)
    return 0.0;
  return path_measurement_contribution_dx(path, 0, path->length-1);
}

static void combine(reservoir_t *s, const reservoir_t *r) {
  // can't combine reservoir with itself (happens when random_neighbor fails for example)
  if(s == r) { printf("tried to combine reservoir with itself\n"); return; } 

  // reservoirs can't be empty (although the path in the reservoir can still be the null sample (= uninitialized path))
  if(s->c <= 0. && r->c <= 0.) { printf("tried to combine empty reservoirs"); return;}

  md_t w_r = 0.;
  md_t w_s = 0.;
  md_t phat_r = p_hat(r->path);
  md_t phat_s = p_hat(s->path);
  //double rm = r->c * phat_r;
  //double sm = s->c * phat_s;

  // bias! how?
  //double mis_r = (rm + sm) > 0. ? rm/(rm+sm) : 0.;
  //double mis_s = (rm + sm) > 0. ? sm/(rm+sm) : 0.;
  //printf("rm %f, sm %f, mis_r %f, mis_s %f\n", rm, sm, mis_r, mis_s);
  double mis_r = r->c / (s->c + r->c);
  double mis_s = s->c / (s->c + r->c);


  // if r has a valid sample (otherwise we call p_hat with an uninitialized_path which can give NaN)
  if(r->w_sum > 0.) {
    // TODO: construct path  [ s->v[0], s->v[1], r->v[2] ]
    w_r = mis_r * phat_r * r->W; // assume s and r for same pixel (starting vertex), then f(r.Y) * r.W = r.w_sum
    //w_r = mis_r * path_throughput(s->v[0], s->v[1], r->v[2]) * r->W; // f_q(r.Y) * estimate for 1/p_q'(r.Y)
  }
  
  // if s has a valid sample
  if(s->w_sum > 0.)
    w_s = mis_s * phat_s * s->W;

  s->w_sum = w_s;
  update(s, r->path, w_r, r->c);

  // update estimator
  if(s->w_sum > 0.) {
    s->W = s->w_sum / p_hat(s->path);
  }

  // confidence capping
  if(s->c > 20.) s->c = 20.;
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
      r->w_sum = 0.;
      r->c = 0.;
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
  assert(!(init_path->v[0].mode & s_emit));

  // reset
  r->c = 0.;
  r->w_sum = 0.;
  r->W = 0.;

  const int M = 8;
  path_t path;
  for(int k = 0; k < M; k++) {
    path_copy(&path, init_path);
    
    // direct illumination, fails when hitpoint of init_path on envmap
    if(nee_sample(&path)) continue;

    md_t w = 0.0;
    md_t f = p_hat(&path);
    md_t p = path.v[2].pdf;
    if(p > 0.0)
      w = (1.0/M) * (f/p); // 1/M uniform weights

    update(r, &path, w, 1.); // new independent sample gets confidence = 1
  }

  // update contribution weight W (= estimator for 1/p(r.Y)), only fails if all M samples were 0 samples
  if(r->w_sum > 0) {
    r->W = r->w_sum / p_hat(r->path);
  }
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

  // combine with existing reservoir
  combine(r, &rris);

  // all 0 samples
  if(r->w_sum <= 0.) {
    free(rris.path);
    return; 
  }

  // lambda shifting
  // lost code...
  
  // estimator f(r.Y) * r.W
    // quick fix: multiply by 1/v[0].pdf * 1/v[1].pdf, because r.W is an estimator of only 1/v[2].pdf!
  pointsampler_splat(r->path, md_2f(f(r->path) * r->W / (path->v[0].pdf * path->v[1].pdf)));
  
  path_copy(path, r->path);

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
