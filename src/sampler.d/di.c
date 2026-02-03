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

// set pixel of paths ourself to control/disable anti-aliasing (necessary because we use only one reservoir per pixel)
static void get_pixel(uint64_t index, uint64_t *i, uint64_t *j) {
  // uint64_t w = view_width();
  // uint64_t h = view_height();
  // *i = (index / h) % w;
  // *j = index % h;
  // halton sequence
  float pixel_i, pixel_j;
  pointsampler_pixel(index, &pixel_i, &pixel_j);
  *i = (uint64_t)(pixel_i);
  *j = (uint64_t)(pixel_j);
}

void sampler_create_path(path_t *path)
{  
  // get pixel from path index
  //uint64_t i, j;
  //get_pixel(path->index, &i, &j);
  
  // extend path once to determine pixel on camera and first vertex
  //path_init(path, path->index, path->sensor.camid);
  //path_set_pixel(path, (float)i+0.5f, (float)j+0.5f); // +0.5 for center of pixel (no anti-aliasing!)

  if(path_extend(path)) return;
  if(path_extend(path)) return;

  // direct illumination
  if(nee_sample(path)) return;

  //if(!path_visible(path, 3)) return;

  pointsampler_splat(path, path_throughput(path));
  return;

  md_t f = path_measurement_contribution_dx(path, 0, path->length-1);
  md_t pdf = path_pdf(path);

  if(pdf <= 0.)
    return;

  pointsampler_splat(path, md_2f(f / pdf));
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
