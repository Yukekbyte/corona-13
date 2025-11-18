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

#include "sampler.h"
#include "spectrum.h"
#include "pointsampler.h"

// resampled importance sampling

typedef struct reservoir_t {
  path_t *path; // output sample
  double w_sum; // sum of weights
  uint64_t M; // number of samples seen
} reservoir_t;

typedef struct sampler_t {} sampler_t;

int update(reservoir_t *r, path_t *path, double weight) {
  r->w_sum += weight;
  r->M++;
  if (((double)rand()/(double)RAND_MAX) < weight / r->w_sum) {
    r->path = path;
    return 1;
  }
  return 0;
}
    
sampler_t *sampler_init() {return 0;}
void sampler_cleanup(sampler_t *s) {}
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
  // keep this stuff on the stack (for now)
  reservoir_t r;
  path_t p;
  r.path = &p;
  r.w_sum = 0;
  r.M = 0;

  // streaming RIS
  int M = 8;
  for(int i = 0; i < M; i++) {
    path_init(&p, path->index, path->sensor.camid);

    // direct illumination
    if(path_extend(&p)) continue;
    if(path_extend(&p)) continue;

    double w = 1.0/M * path_throughput(&p);
    update(&r, &p, w);
  }

  // call pointsampler_splat() on chosen sample
  pointsampler_splat(r.path, path_throughput(r.path) * r.w_sum);

  // copy path at the end
  path_copy(path, r.path);
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
