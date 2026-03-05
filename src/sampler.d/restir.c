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

#define M 8
#define NEIGHBOUR_COUNT 4
#define NEIGHBOUR_RADIUS 10 // radius must be sufficiently big for the neighbour count.
#define PAIRWISE_COMBINE 1
#define SPATIAL_REUSE_PASSES 1
#define TEMPORAL_REUSE 0
#define CONFIDENCE_CAP 100. // is a double

// Reservoir-based Spatio-Temporal Importance Resampling (ReSTIR)

typedef struct reservoir_t {
  path_t *path; // output sample
  double w_sum; // sum of weights
  double c;     // confidence weight of output (= the amount of samples behind the output sample)
  double W;     // contribution weight: estimate for 1/p
} reservoir_t;

typedef struct pixel_t {
  uint64_t i;
  uint64_t j;
  float _i;
  float _j;
} pixel_t;

typedef struct sampler_t {
  reservoir_t **reservoirs;
  int spatial_reuse_passes;
} sampler_t;

// returns random double between 0 and 1
// separate from pointsampler for things that should always be random (never stratified)
static double random_uniform() {
  return ((double)rand()/(double)RAND_MAX);
}

float randf() { 
  return rand() / (float)RAND_MAX; 
}

// fractional part of float
static inline float fractf(float x) { return x - floorf(x); }

static void get_pixel(const uint64_t index, pixel_t *q) {
  pointsampler_pixel(index, &q->i, &q->j, &q->_i, &q->_j);
}

static void get_pixel_linear(const uint64_t index, pixel_t *q) {
  pointsampler_pixel_linear(index, &q->i, &q->j, &q->_i, &q->_j);
}

// Updates reservoir with sample and weight.
static void update(reservoir_t *r, path_t *path, double weight, double c) {
  r->w_sum += weight;
  r->c += c;
  if (random_uniform() * r->w_sum < weight)
    path_copy(r->path, path);
}

// Use the integrand f as target function p_hat
static double p_hat(path_t *path) {
  if(is_null(path))
    return 0.0;
  return md_hsum(path_measurement_contribution_dx(path, 0, path->length-1));
}

// Perform Resampled Importance Sampling (streaming RIS)
// r will be reset and filled with initial samples for pixel q
static void ris(pixel_t q, reservoir_t *r) {
  // reset
  set_null(r->path);
  r->c = 0.;
  r->w_sum = 0.;
  r->W = 0.;

  for(int i = 0; i < M; i++) {
    path_t path;
    path_init(&path, 0, 0);
    path_set_pixel(&path, q._i, q._j);
    
    if(path_extend(&path)) break;

    if(path.v[path.length-1].flags & s_environment) {
      // give reservoir non-zero confidence to serve as indicator
      r->c = 1;
      path_copy(r->path, &path);
      break;
    }
    
    // generate path tree
    while(1) {
      // sample light source
      if(nee_sample(&path)) break; // breaks when envmap is hit or path becomes too long

      double phat = p_hat(&path);
      double pdf = md(path_pdf(&path), 0); // hero wavelength pdf
      if(phat > 0. && pdf > 0.)
        update(r, &path, phat/pdf, 1.); // mis weight is 1 for samples of same path tree because each sample has a different length
      else
        r->c += 1.;
      path_pop(&path);
      
      // extend path
      if(path_extend(&path)) break;  
    }
  }

  // mis weights for each chosen sample per path tree is 1/M (
  // added now instead of before for simplicity, is mathematically equivalent
  r->w_sum *= 1./M;

  // update contribution weight W (= estimator for 1/p(r.Y)), only fails if all M samples were 0 samples
  if(not_null(r->path)) {
    r->W = r->w_sum / p_hat(r->path);
  } 
}

float shift(path_t *shifted, pixel_t q, const path_t *source_path) {
  if(is_null(source_path)) {
    set_null(shifted);
    return 0.0;
  }

  float J = path_shift(shifted, q._i, q._j, source_path, 1);
  
  // check if shift failed
  if (J == 0.0f || p_hat(shifted) == 0.0f) {
    set_null(shifted);
    return 0.0;
  }

  return J;
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

  // MIS weights
  if(is_null(s->path))
    mis_s = 0.0f;
  else
    mis_s = s->c * p_hat(s->path) / (s->c * p_hat(s->path) + r->c * p_hat(&s_path_from_qr) * Js);
  
  if(is_null(r->path))
    mis_r = 0.0f;
  else
    mis_r = r->c * p_hat(r->path) / (r->c * p_hat(r->path) + s->c * p_hat(&r_path_from_qs) * Jr);
  
  // resampling weights
  double w_r = mis_r * p_hat(&r_path_from_qs) * r->W * Jr;
  double w_s = mis_s * p_hat(s->path) * s->W;

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
double p_hat_from(path_t *y, pixel_t q) {
  path_t x;
  float J = shift(&x, q, y);
  if(not_null(&x))
    return p_hat(&x) * J;
  return 0;
}

double p_hat_from_opt(path_t *x, float J) {
  return p_hat(x) / J;
}

double mis(path_t *x, path_t *y, float J, double cy, pixel_t q[], double c[]) {
  if(is_null(y)) return 0.0;
  
  double num = cy * p_hat_from_opt(x, J);
  if(num <= 0.0) return 0.0;
  
  double denom = num;
  for(int i = 0; i < NEIGHBOUR_COUNT; i++)
    denom += c[i] * p_hat_from(y, q[i]);

  return num / denom;
}

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
// Combine reservoir s with r.
// Assumes s and r come from the same domain (pixel)
static void combine_temporal(reservoir_t *s, const reservoir_t *r) {
  // don't combine reservoir with itself (happens when random_neighbor fails and returns its own reservoir)
  if(s == r) { printf("tried to combine reservoir with itself\n"); return; } 

  // reservoirs can't be empty (although the path in the reservoir can still be a null sample)
  if(s->c <= 0. && r->c <= 0.) { printf("tried to combine empty reservoirs\n"); return; }
  
  // Shift lambda of r, should be deterministic
  
  r->path->lambda = spectrum_sample_lambda(pointsampler(r->path, s_dim_lambda), NULL);

  // MIS weights
  double total = (s->c + r->c);
  double mis_s = s->c / total;
  double mis_r = r->c / total;
  
  // resampling weights
  md_t w_r = mis_r * p_hat(r->path) * r->W;
  md_t w_s = mis_s * p_hat(s->path) * s->W;

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

// Pick neighbours with R2 sequence
// Partly ChatGPT
static void random_neighbors(pixel_t q, const path_t *path, reservoir_t **ns, pixel_t *qns)
{
    const uint64_t w = view_width();
    const uint64_t h = view_height();

    const int d = NEIGHBOUR_RADIUS;
    const int side = 2 * d + 1;
    const int area = side * side;

    // R2 constants
    const float a1 = 0.7548776662466927f;  // 1/phi
    const float a2 = 0.5698402909980532f;  // 1/phi^2
    float seed_u = randf();
    float seed_v = randf();

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
  new.path = (path_t *)malloc(sizeof(path_t));
  // Intial RIS
  ris(q, &new);
  // combine with temporal neighbour (previous reservoir)
  combine_temporal(r, &new);
  free(new.path);
  #else
  // Intial RIS
  ris(q, r);
  #endif
}

int sampler_passes() { return SPATIAL_REUSE_PASSES; }
void sampler_pass_sample(uint64_t index) {
  // get pixels with pointsampler sequence to avoid artifacts and race conditions
  pixel_t q;
  get_pixel(index, &q);
  reservoir_t *r = &rt.sampler->reservoirs[q.i][q.j];

  // check for env map hit
  if(not_null(r->path) && r->path->v[r->path->length-1].flags & s_environment) return;

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

void sampler_create_path(path_t *path)
{
  pixel_t q;
  get_pixel_linear(path->index, &q);
  reservoir_t *r = &rt.sampler->reservoirs[q.i][q.j];

  // don't splat null sample
  if(is_null(r->path)) return;

  // estimator f(r.Y) * r.W
  const md_t estimator = md_mul(path_measurement_contribution_dx(r->path, 0, r->path->length-1), md_set1(r->W));
  const mf_t w = path_pdf_hero(r->path);

  pointsampler_splat(r->path, mf_mul(w, md_2f(estimator)));
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : restir\n");
}
