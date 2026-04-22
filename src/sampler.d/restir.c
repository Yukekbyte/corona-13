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
#include "reservoir.h"
#include "gris.h"

#define SPATIAL_REUSE_PASSES 3
#define SPLAT_COUNT_CORRECTOR mf_set1(1.0/((double)(M + SPATIAL_REUSE_PASSES + 1)))

// Reservoir-based Spatio-Temporal Importance Resampling (ReSTIR)

typedef struct sampler_t {
  reservoir_t **reservoirs;
  int spatial_reuse_passes;
} sampler_t;

// fractional part of float
static inline float fractf(float x) { return x - (int)x; }

static void get_pixel(const uint64_t index, pixel_t *q) {
  pointsampler_pixel(index, &q->i, &q->j, &q->_i, &q->_j);
}

static void get_pixel_linear(const uint64_t index, pixel_t *q) {
  pointsampler_pixel_linear(index, &q->i, &q->j, &q->_i, &q->_j);
}

static inline void splat(path_t *p, mf_t estimator) {
  if(mf_any(mf_gt(estimator, mf_set1(0.0)))) {
    pointsampler_splat(p, mf_mul(estimator, SPLAT_COUNT_CORRECTOR));
  }
}

// Pick neighbours with R2 sequence
// Partly ChatGPT
static void random_neighbors(pixel_t q, const path_t *path, reservoir_t **ns, pixel_t *qns)
{
  const int tid = common_get_threadid();

  const uint64_t w = view_width();
  const uint64_t h = view_height();

  const int d = NEIGHBOUR_RADIUS;
  const int side = 2 * d + 1;
  const int area = side * side;

  // R2 constants
  const float a1 = 0.7548776662466927f;  // 1/phi
  const float a2 = 0.5698402909980532f;  // 1/phi^2
  float seed_u = points_rand(rt.points, tid);
  float seed_v = points_rand(rt.points, tid);

  int found = 0;
  int k     = 0;
  const int MAX_ATTEMPTS = area * 2;

  while(found < NEIGHBOUR_COUNT && k < MAX_ATTEMPTS)
  {
    // 2D R2 sequence in [0,1)^2
    float u = fractf((k + 1) * a1 + seed_u);
    float v = fractf((k + 1) * a2 + seed_v);

    // Map to square neighborhood
    int dx = (int)(u * side) - d;
    int dy = (int)(v * side) - d;

    int l = (int)q.i + dx;
    int m = (int)q.j + dy;

    // Clamp to image bounds
    l = CLAMP(l, 0, (int)w - 1);
    m = CLAMP(m, 0, (int)h - 1);

    // Skip self
    if((uint64_t)l == q.i && (uint64_t)m == q.j) {
      k++;
      continue;
    }

    // Skip duplicates
    int duplicate = 0;
    for(int i = 0; i < found; i++) {
      if(qns[i].i == (uint64_t)l &&
      qns[i].j == (uint64_t)m) {
        duplicate = 1;
        break;
      }
    }

    if(duplicate) {
      k++;
      continue;
    }

    pixel_t *qn = &qns[found];
    qn->i = (uint64_t)l;
    qn->j = (uint64_t)m;
    reservoir_t *n = &rt.sampler->reservoirs[qn->i][qn->j];
    
    // Check for geometric similarity (up to a point)
    if(k <= 3*NEIGHBOUR_COUNT && k <= 0.5*MAX_ATTEMPTS) {
      // angle between normals < 25 deg (0.435 rad) 
      if(acosf(dotproduct(path->v[1].hit.n, n->path->v[1].hit.n)) > 0.436f) {
        k++;
        continue;
      }
    
      // depth difference can't be more than 10 percent
      float depthratio = path->e[1].dist / n->path->e[1].dist;
      if(depthratio < 0.9f || 1.1f < depthratio) {
        k++;
        continue;
      }
    }

    pointsampler_subpixel(qn->i, qn->j, &qn->_i, &qn->_j);

    ns[found] = n;
    
    found++;
    k++;
  }

  // Sanity check
  assert(found == NEIGHBOUR_COUNT);
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
      r->envmap = 0;
    }
  }

  s->spatial_reuse_passes = SPATIAL_REUSE_PASSES;

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

void sampler_prepare_sample(uint64_t index) {
  // get pixels linearly
  pixel_t q;
  get_pixel_linear(index, &q);
  reservoir_t *r = &rt.sampler->reservoirs[q.i][q.j];
  
  #if TEMPORAL_REUSE
    reservoir_t new;
    path_t path;
    new.path = &path;
    // Intial RIS
    ris(q, &new, NULL);

    if(new.envmap)
      r->envmap = 1; 
    else // combine with temporal neighbour (previous reservoir)
      combine_temporal(r, &new);

  #else
    // Intial RIS
    ris(q, r, NULL);
  #endif
}

int sampler_passes() { return SPATIAL_REUSE_PASSES; }
void sampler_pass_sample(uint64_t index) {
  // get pixels with pointsampler sequence to avoid artifacts and race conditions
  pixel_t q;
  get_pixel(index, &q);
  reservoir_t *r = &rt.sampler->reservoirs[q.i][q.j];

  if(r->envmap) return;

  // Spatial re-use pass
  reservoir_t *ns[NEIGHBOUR_COUNT];
  pixel_t qns[NEIGHBOUR_COUNT];
  random_neighbors(q, r->path, ns, qns);
  
  #if PAIRWISE_COMBINE
    for(int k = 0; k < NEIGHBOUR_COUNT; k++)
      combine_pair(q, r, qns[k], ns[k]);
  #else
    combine(q, r, qns, ns);
  #endif

  // splat sample anyway
  // if(is_null(r->path) || r->envmap) return;
  // mf_t estimator = md_2f(md_mul(path_measurement_contribution_dx(r->path, 0, r->path->length-1), md_set1(r->W)));
  // splat(r->path, estimator);
}

void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

void sampler_create_path(path_t *path)
{
  pixel_t q;
  get_pixel_linear(path->index, &q);
  reservoir_t *r = &rt.sampler->reservoirs[q.i][q.j];

  // don't splat null sample and don't support envmaps
  if(is_null(r->path) || r->envmap) return;

  // estimator f(r.Y) * r.W
  const mf_t estimator = md_2f(md_mul(path_measurement_contribution_dx(r->path, 0, r->path->length-1), md_set1(r->W)));
  if(mf_any(mf_gt(estimator, mf_set1(0.0)))) {
    pointsampler_splat(r->path, estimator);
  }
  // splat(r->path, estimator);
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : restir\n");
}
