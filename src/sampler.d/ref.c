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


typedef struct sampler_t {} sampler_t;

sampler_t *sampler_init() {return 0;}
void sampler_cleanup(sampler_t *s) {}
void sampler_prepare_frame(sampler_t *s) {}
void sampler_clear(sampler_t *s) {}

static void get_pixel_linear(const uint64_t index, pixel_t *q) {
  pointsampler_pixel_linear(index, &q->i, &q->j, &q->_i, &q->_j);
}

void sampler_create_path(path_t *path)
{
  pixel_t q;

  get_pixel_linear(path->index, &q);
  path_init(path, path->index, 0);
  path_set_pixel(path, q._i, q._j);

  // path->index = (uint32_t)(points_rand(rt.points, common_get_threadid()) * 4294967296.0f); // 2^32

  // extend path once to determine pixel on camera and first vertex (hitpoint)
  if(path_extend(path)) return;

  while(1) {
    // sample light source
    if(nee_sample(path)) break;

    md_t measurement = path_measurement_contribution_dx(path, 0, path->length-1);
    double pdf = md_hsum(path_pdf(path));
    pointsampler_splat(path, md_2f(md_div(measurement, md_set1(pdf))));

    // pop area sampled light vertex
    path_pop(path);

    // extend path with bsdf sampling
    if(path_extend(path)) break;
  }
}

mf_t sampler_throughput(path_t *path)
{
  if(path->length < 1) return mf_set1(0.0f);
  const md_t measurement = path_measurement_contribution_dx(path, 0, path->length-1);
  if(mf_all(mf_lte(md_2f(measurement), mf_set1(0.0f))))
    return mf_set1(0.0f);
  md_t pdf = md_set1(1.0);
  for(int k=0;k<path->length;k++)
    pdf = md_mul(pdf, mf_2d(path_pdf_extend(path, k)));
  return md_2f(md_div(measurement, pdf));
}

void sampler_print_info(FILE *fd)
{
  fprintf(fd, "sampler  : ref\n");
}
