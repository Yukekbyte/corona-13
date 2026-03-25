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

// resampled importance sampling

typedef struct sampler_t {}
sampler_t;

static void get_pixel_linear(const uint64_t index, pixel_t *q) {
  pointsampler_pixel_linear(index, &q->i, &q->j, &q->_i, &q->_j);
}

sampler_t *sampler_init() { return 0; }

void sampler_cleanup(sampler_t *s) {}
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

static inline void splat(path_t *p, mf_t estimator) {
  // if(mf_any(mf_gt(estimator, mf_set1(0.0)))) {
  //   const mf_t w = path_pdf_hero(p);
  //   //pointsampler_splat(p, mf_mul(w, estimator));
  //   pointsampler_splat(p, mf_mul(w, mf_mul(estimator, mf_set1(1/(M+1)))));
  // }
}

void sampler_create_path(path_t *path)
{  
  // get pixel & reservoir from path index
  reservoir_t r;
  pixel_t q;
  path_t p;
  r.path = &p;

  get_pixel_linear(path->index, &q);
  
  // // inital candidate generation
  ris(q, &r, splat);

  // don't splat null sample
  if(is_null(r.path) || r.envmap) return;

  // estimator f(r.Y) * r.W
  const mf_t estimator = md_2f(md_mul(path_measurement_contribution_dx(r.path, 0, r.path->length-1), md_set1(r.W)));
  if(mf_any(mf_gt(estimator, mf_set1(0.0)))) {
    const mf_t w = path_pdf_hero(r.path);
    //pointsampler_splat(p, mf_mul(w, estimator));
    pointsampler_splat(r.path, mf_mul(w, estimator));
  }
  // const mf_t estimator = md_2f(md_mul(path_measurement_contribution_dx(r.path, 0, r.path->length-1), md_set1(r.W)));
  // splat(r.path, estimator);
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : ris\n");
}
