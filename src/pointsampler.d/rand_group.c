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
#include "view.h"
#include "assert.h"

#define S 100 // Distance between pixels when sampling the sensor

typedef struct fake_randoms_t
{
  int enabled;
  float rand[40];
}
fake_randoms_t;

typedef struct group_t
{
  uint64_t size;
  uint64_t x, y;
  uint64_t count_x, count_y;
} 
group_t;

typedef struct pointsampler_t
{
  fake_randoms_t *rand;
  group_t group[S*S];
  uint64_t cum_group_size[S*S];
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

  uint64_t w = view_width();
  uint64_t h = view_height();

  uint64_t i = 0;
  group_t *g;
  for(uint64_t x = 0; x < S; x++) 
    for(uint64_t y = 0; y < S; y++){
      g = &s->group[i];
      assert(x < w && y < h);
      g->count_x = (w - 1 - x) / S + 1;
      g->count_y = (h - 1 - y) / S + 1;
      g->size = g->count_x * g->count_y;
      g->x = x;
      g->y = y;
      s->cum_group_size[i] = g->size + (i > 0 ? s->cum_group_size[i-1] : 0);
      i++;
  }

  return s;
}

float pointsampler(path_t *p, int dim)
{
  // pure random mersenne twister
  const int tid = common_get_threadid();
  if(rt.pointsampler->rand[tid].enabled)
    return rt.pointsampler->rand[tid].rand[dim];
  return points_rand(rt.points, tid);
}

static inline uint64_t find_group(const pointsampler_t *s, uint64_t index)
{
  uint64_t lo = 0;
  uint64_t hi = S*S-1; // last valid group index

  while (lo < hi)
  {
    uint64_t mid = (lo + hi) >> 1;

    if (index < s->cum_group_size[mid])
      hi = mid;
    else
      lo = mid + 1;
  }

  return lo;
}

void pointsampler_pixel_linear(uint64_t index, uint64_t *x, uint64_t *y, float *pixel_i, float *pixel_j)
{
  uint64_t width  = view_width();
  uint64_t height = view_height();
  uint64_t total = width * height;

  uint64_t wrapped = index % total;

  *x = wrapped % width;
  *y = wrapped / width;

  pointsampler_subpixel(*x, *y, pixel_i, pixel_j);
}

void pointsampler_pixel(uint64_t index, uint64_t *x, uint64_t *y, float *pixel_i, float *pixel_j)
{
  const pointsampler_t *s = rt.pointsampler;

  uint64_t wrapped = index % s->cum_group_size[S*S-1];
  uint64_t g = find_group(s, wrapped);
  
  uint64_t group_start = (g == 0) ? 0 : s->cum_group_size[g - 1];

  uint64_t local = wrapped - group_start;
  
  const group_t *grp = &s->group[g];

  uint64_t tx = local % grp->count_x;
  uint64_t ty = local / grp->count_x;

  *x = grp->x + tx * S;
  *y = grp->y + ty * S;

  pointsampler_subpixel(*x, *y, pixel_i, pixel_j);
}

void pointsampler_subpixel(uint64_t x, uint64_t y, float *pixel_i, float*pixel_j) {
  *pixel_i = (float)x + 0.5f;
  *pixel_j = (float)y + 0.5f;
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
