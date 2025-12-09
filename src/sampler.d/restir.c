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
#include "pathspace/multichain.h"
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
  double c;     // confidence weight of output (= the amount of samples behind the output sample)
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

// Combine reservoir r with s.
// Path is a length 2 path that contains the camera vertex and hitpoint for which the throughput in the reservoirs will be calculated.
// Reservoir s is assumed to already be based on the given path.
static void combine(const path_t *path, reservoir_t *s, const reservoir_t *r, const uint64_t i, const uint64_t j) {
  // can't combine reservoir with itself (happens when random_neighbor fails for example)
  if(s == r) { printf("tried to combine reservoir with itself\n"); return; } 

  // reservoirs can't be empty (although the path in the reservoirs can still be the null sample)
  if(s->c <= 0. && r->c <= 0.) { printf("tried to combine empty reservoirs"); return;}

  // shift r->path
  path_t shifted;
  path_copy(&shifted, path);

  if(not_null(r->path)) {
    // if(i == 800 && j == 200) {
    //   int status = path_shift(&shifted, r->path); // s.v[0] -> s.v[1] -> r.v[2]
    //   printf("shift status: %d | p_hat(&shifted): %f\n", status, p_hat(&shifted));
    // }
    // path_copy(&shifted, r->path); 
    multichain_perturb_connect(r->path, &shifted, 1);
    // transform probability to on-surface probability at vertex v
    shifted.v[1].pdf = mf_mul(shifted.v[1].pdf, path_G(&shifted, 1));
    if(i == 800 && j == 200){
      printf("------(1) shifted->v[1].pdf = %f\n", shifted.v[1].pdf);
    }

    if(i == 800 && j == 200){
      printf("--- r->path ---\n");
      path_print(r->path, stdout);
      // printf("--- s->path ---\n");
      // path_print(s->path, stdout);
      printf("--- shifted ---\n");
      path_print(&shifted, stdout);
      printf("\n");

      printf("r->path Estimator: %f, f: %f, W: %f, 1/v[0].pdf: %f, 1/v[1].pdf: %f\n", md_2f(f(r->path) * r->W / (r->path->v[0].pdf * r->path->v[1].pdf)),
                                                                              f(r->path), r->W, 1/r->path->v[0].pdf, 1/r->path->v[1].pdf);
    }
  }
  else {
    if(i == 800 && j == 200){
      printf("r->path is null\n");
    }
    set_null(&shifted);
  }
  
  // calculate weights
  double mis_s = s->c / (s->c + r->c);
  double mis_r = r->c / (s->c + r->c);
  md_t w_r = mis_r * p_hat(&shifted) * r->W;
  md_t w_s = mis_s * p_hat(s->path) * s->W;
  
  // combine reservoirs in s
  s->w_sum = w_s;
  int stat = update(s, &shifted, w_r, r->c);

  // update estimator
  if(not_null(s->path))
    s->W = s->w_sum / p_hat(s->path); // NaN if s->path is null (not really a problem)
  else
    s->W = 0.;

  if(i == 800 && j == 200) {

    if(stat) {
      printf("picked shifted\n");
    } else {
      printf("picked prev\n");
    }
    printf("final s Estimator: %f, f: %f, W: %f, 1/v[0].pdf: %f, 1/v[1].pdf: %f\n\n", md_2f(f(s->path) * s->W / (s->path->v[0].pdf * s->path->v[1].pdf)),
    f(s->path), s->W, 1/s->path->v[0].pdf, 1/s->path->v[1].pdf);
  }
  
  // confidence capping
  if(s->c > 20.) s->c = 20.;
}

static reservoir_t *random_neighbor(uint64_t i, uint64_t j, path_t *path) {
  uint64_t w = view_width();
  uint64_t h = view_height();

  int l, m; // int instead of uint64_t so they can be negative (important for clamping)
  reservoir_t *r;

  const int d = 30; // sample in 20x20 square around (i, j)
  const int MAX_ATTEMPTS = 10; // when i == l && j == m or neighbor not geometrically similar enough
  for(int k = 0; k < MAX_ATTEMPTS; k++) {
    // TODO: make edge clamping uniform probablility
    l = i + (2*d*random_uniform() - d);
    m = j + (2*d*random_uniform() - d);

    l = CLAMP(l, 0, w-1);
    m = CLAMP(m, 0, h-1);

    uint64_t ni = (uint64_t)l;
    uint64_t nj = (uint64_t)m;

    assert(ni < w && nj < h);

    if(ni == i && nj == j) continue;

    r = &rt.sampler->reservoirs[ni][nj];
    
    // geometric similarity
    // TODO: r->path->v[1] must be similar to path->v[1]
    return r;
  }

  return &rt.sampler->reservoirs[i][j]; // return same reservoir
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
  assert(!(init_path->v[0].mode & s_emit));

  // reset
  set_null(r->path);
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
  
  // extend path once for camera vertex and first hitpoint
  path_init(path, path->index, path->sensor.camid);
  path_set_pixel(path, (float)i+0.5f, (float)j+0.5f); // +0.5 for center of pixel (no anti-aliasing!)
  if(path_extend(path)) return;

  // check for env map hit
  if(path->v[path->length-1].flags & s_environment) {
    pointsampler_splat(path, path_throughput(path));
    return;
  }

  // inital candidate generation
  reservoir_t rris;
  rris.path = (path_t*)malloc(sizeof(path_t));
  ris(&rris, path);

  // combine with existing reservoir
  combine(path, r, &rris, i, j);

  // spatial re-use
  const int neighbors = 3;
  for(int k = 0; k < neighbors; k++) {
    // TODO: pick neighbors without replacement
    if(i == 800 && j == 200)
      printf("+++ Combining with neighbor +++\n");
    combine(path, r, random_neighbor(i, j, r->path), i, j);
  }

  // don't splat null sample
  if(is_null(r->path)) {
    if(i == 800 && j == 200) printf("Estimator: / (r->path is null)\n\n");
    free(rris.path);
    return; 
  }

  // estimator f(r.Y) * r.W
    // multiply by 1/v[0].pdf * 1/v[1].pdf, because r.W is an estimator of only 1/v[2].pdf!
  if(i == 800 && j == 200) {
    printf("Estimator: %f, f: %f, W: %f, 1/v[0].pdf: %f, 1/v[1].pdf: %f\n\n", md_2f(f(r->path) * r->W / (r->path->v[0].pdf * r->path->v[1].pdf)),
                                                                              f(r->path), r->W, 1/r->path->v[0].pdf, 1/r->path->v[1].pdf);
  }
  pointsampler_splat(r->path, md_2f(f(r->path) * r->W / (r->path->v[0].pdf * r->path->v[1].pdf)));
  
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
