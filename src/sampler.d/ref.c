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

static void get_pixel_linear(uint64_t index, uint64_t *i, uint64_t *j, float *_i, float *_j) {
  pointsampler_pixel_linear(index, i, j, _i, _j);
}

void sampler_create_path(path_t *path)
{  
  // get pixel from path index
  uint64_t i, j;
  float _i, _j;
  get_pixel_linear(path->index, &i, &j, &_i, &_j);
  path_set_pixel(path, _i, _j);
  
  // extend path once to determine pixel on camera and first vertex
  if(path_extend(path)) return;

  const int max_length = 10;
  while(path->length < max_length) {
    // sample light source
    if(nee_sample(path)) break;

    // path.throughput == p_hat/pdf
    pointsampler_splat(path, path->throughput);
    path_pop(path);

    // extend path
    if(path_extend(path)) break;
  }
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
  fprintf(fd, "sampler  : ref\n");
}
