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
    along with corona-13. If not, see <http://www.gnu.org/licenses/>.
*/

#include "pointsampler.h"
#include "sampler.h"
#include "render.h"
#include "points.h"
#include "threads.h"
#include "ext/halton/halton.h"

typedef struct fake_randoms_t
{
  int enabled;
  float rand[40];
}
fake_randoms_t;

typedef struct pointsampler_t
{
  fake_randoms_t *rand;
}
pointsampler_t;

void pointsampler_print_info(FILE *f)
{
  fprintf(f, "mutations: none\n");
}

pointsampler_t *pointsampler_init(uint64_t frame)
{
  pointsampler_t *s = calloc(1, sizeof(*s));
  s->rand = calloc(rt.num_threads, sizeof(*s->rand));
  return s;
}

float pointsampler(path_t *p, int dim)
{
  const int tid = common_get_threadid();
  if(rt.pointsampler->rand[tid].enabled)
    return rt.pointsampler->rand[tid].rand[dim];
  
  // store wavelength random number because we want equally spaced wavelengths.
  // The calls to pointsampler with dim_lambda on the same path in pathspace extend must give the same random number.
  if(dim == s_dim_lambda) {
    if(p->lambda_rand >= 0.0f)
      return p->lambda_rand;
    else {
      float r = points_rand(rt.points, tid);
      p->lambda_rand = r; 
      return r;
    }
  }
  
  // pure random mersenne twister
  return points_rand(rt.points, tid);
}

// A variant of pointsampler that stores the random numbers.
// Used for 3 dimensions during the random replay shifts in pathspace.
float pointsampler_store(path_t *p, int v, int dim)
{
  // pure random mersenne twister
  const int tid = common_get_threadid();

  float r;
  if(rt.pointsampler->rand[tid].enabled)
    r = rt.pointsampler->rand[tid].rand[dim];
  else 
    r = points_rand(rt.points, tid);
  
  if(dim == s_dim_omega_x)
    p->v[v].rands.omega_x = r;
  else if(dim == s_dim_omega_y)
    p->v[v].rands.omega_y = r;
  else if(dim == s_dim_scatter_mode)
    p->v[v].rands.scatter_mode = r;
  else
    printf("Random number storage not supported for this dimension\n");

  return r;
}

// Gets the pixels linearly.
void pointsampler_pixel_linear(uint64_t index, int *x, int *y)
{
  uint64_t width  = view_width();
  uint64_t height = view_height();
  uint64_t total = width * height;

  uint64_t wrapped = index % total;

  *x = (int)(wrapped % width);
  *y = (int)(wrapped / width);
}

// Randomize within one pixel for anti-aliasing.
void pointsampler_subpixel(int x, int y, float *pixel_i, float*pixel_j) {
  // random uniform within pixel
  const int tid = common_get_threadid();
  *pixel_i = (float)x + points_rand(rt.points, tid);
  *pixel_j = (float)y + points_rand(rt.points, tid);
}

void pointsampler_splat(path_t *p, mf_t value)
{
  render_splat(p, value);
}

void pointsampler_mutate(path_t *curr, path_t *tent)
{
  path_init(tent, tent->index, tent->sensor.camid);
  sampler_create_path(tent);
}

void pointsampler_mutate_with_pixel(path_t *curr, path_t *tent, float i, float j)
{
  path_init(tent, tent->index, tent->sensor.camid);
  path_set_pixel(tent, i, j);
  sampler_create_path(tent);
}

int pointsampler_accept(path_t *curr, path_t *tent) { return 0; }
void pointsampler_clear() {}
void pointsampler_cleanup(pointsampler_t *s)
{
  free(s->rand);
  free(s);
}
void pointsampler_set_large_step(pointsampler_t *t, float p_large_step) {}
void pointsampler_reset_thread(pointsampler_t *t) {}
void pointsampler_finalize(pointsampler_t *s) {}

void pointsampler_enable_fake_random(pointsampler_t *s)
{
  const int tid = common_get_threadid();
  s->rand[tid].enabled = 1;
}

void pointsampler_disable_fake_random(pointsampler_t *s)
{
  const int tid = common_get_threadid();
  s->rand[tid].enabled = 0;
}

void pointsampler_set_fake_random(pointsampler_t *s, int dim, float rand)
{
  const int tid = common_get_threadid();
  s->rand[tid].rand[dim] = rand;
}

void pointsampler_prepare_frame(pointsampler_t *s) {}

void pointsampler_stop_learning(pointsampler_t *s) {}
