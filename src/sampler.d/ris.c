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

#include "view.h"
#include "sampler.h"
#include "spectrum.h"
#include "pointsampler.h"

// resampled importance sampling

typedef struct reservoir_t {
  path_t *path; // output sample
  double w_sum; // sum of weights
  uint64_t M; // number of samples seen
} reservoir_t;

typedef struct sampler_t {
  reservoir_t **reservoirs;
} sampler_t;

int update(reservoir_t *r, path_t *path, double weight) {
  r->w_sum += weight;
  r->M++;
  if (((double)rand()/(double)RAND_MAX) < weight / r->w_sum) {
    path_copy(r->path, path);
    return 1;
  }
  return 0;
}
    
sampler_t *sampler_init() {
  int i, j;
  int w = view_width();
  int h = view_height();

  sampler_t *s = (sampler_t *)malloc(sizeof(sampler_t));

  s->reservoirs = (reservoir_t**)malloc(w * sizeof(reservoir_t*));
  for(i = 0; i < w; i++) {
    s->reservoirs[i] = (reservoir_t*)malloc(h * sizeof(reservoir_t));
  }

  reservoir_t *r;
  for(i = 0; i < w; i++) {
    for(j = 0; j < h; j++) {
      r = &s->reservoirs[i][j];
      r->path = NULL;
      r->w_sum = 0;
      r->M = 0;
    }
  }

  return s;
}
void sampler_cleanup(sampler_t *s) {
  int w = view_width();
  
  for(int i = 0; i < w; i++) {
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

void sampler_create_path(path_t *path)
{
  // keep this stuff on the stack
  reservoir_t *r;
  path_t p_init;
  
  // init and extend path once to determine pixel on camera and first vertex (stuff handled by pointsampler)
  path_init(&p_init, path->index, path->sensor.camid);
  if(path_extend(&p_init)) return;

  int i = (int)p_init.sensor.pixel_i;
  int j = (int)p_init.sensor.pixel_j;

  // reset the reservoir (for now)
  r = &rt.sampler->reservoirs[i][j];
  r->path = &p_init;
  r->w_sum = 0;
  r->M = 0;

  // streaming RIS
  int M = 8;
  for(int k = 0; k < M; k++) {
    path_t p; // declaring out of loop creates issues...
    path_copy(&p, &p_init);

    // direct illumination
    // sampling bsdf instead of light source: Fix this!
    if(path_extend(&p)) continue;

    double w = 1.0/M * path_throughput(&p);
    update(r, &p, w);
  }

  // call pointsampler_splat() on chosen sample
  pointsampler_splat(r->path, path_throughput(r->path) * r->w_sum);

  // copy path at the end
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
