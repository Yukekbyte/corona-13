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


typedef struct sampler_t {} sampler_t;

sampler_t *sampler_init() {return 0;}
void sampler_cleanup(sampler_t *s) {}
void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

static inline mf_t sampler_mis(const path_t *p)
{
  md_t pdf = md_set1(1.0);
  for(int v=1;v<p->length;v++)
    pdf = md_mul(pdf, mf_2d(p->v[v].pdf));

  return mf_div(md_2f(pdf), mf_set1(mf_hsum(md_2f(pdf))));
}

static inline mf_t new_lambda(path_t *path) {
  float lf[MF_COUNT];
  mf_t lambda_pdf;
    for(int l=0;l<MF_COUNT;l++)
      lf[l] = fmodf(pointsampler(path, s_dim_lambda) + l/(float)MF_COUNT, 1.0f);
  return spectrum_sample_lambda(mf_loadu(lf), &lambda_pdf);
}

void sampler_create_path(path_t *path)
{  
  // extend path once to determine pixel on camera and first vertex
  if(path_extend(path)) return;

  while(1) {
    // sample light source
    if(nee_sample(path)) break;

    // path.throughput == p_hat/pdf
    if(mf_any(mf_gt(path->throughput, mf_set1(0.0)))) {
      const mf_t w = sampler_mis(path);
      pointsampler_splat(path, mf_mul(w, path->throughput));
    }
    path_pop(path);

    // extend path
    if(path_extend(path)) break;
  }
}

md_t sampler_sum_pdf_dwp(path_t *p)
{
  md_t pdf = md_set1(1.0);
  for(int v=1;v<p->length;v++)
    pdf = md_mul(pdf, mf_2d(mf_div(path_pdf_extend(p, v), mf_set1(path_G(p, v)))));
  return pdf;
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : ref\n");
}
