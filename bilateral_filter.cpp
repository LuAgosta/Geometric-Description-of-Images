/**
 * Copyright (C) 2016, Roberto P.Palomares <roberto.palomares@upf.edu>
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD
 * License. You should have received a copy of this license along
 * this program. If not, see
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */
#include "bilateral_filter.h"


Image<float> BilateralFilter::calculate(FixedImage<float> source,
                                        int r,
                                        int sigma_x,
                                        int sigma_y,
                                        int sigma_c,
                                        int iter)
{
    //TODO: Sin control de parametros
  Image<float> output = source;
  int channels = source.get_number_of_channels();
  int w = source.get_size_x();
  int h = source.get_size_y(); 

  // printf("C:%d W:%d H:%d \n",channels, w, h);


  for (int k = 0; k < iter; k++){
    // printf("Iter: %d\n", k);
    for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      // printf("(%d,%d)\n", x,y);
      int x_min = ((x - r) >= 0) ? (x - r): 0;
      int y_min = ((y - r) >= 0) ? (y - r): 0;
      int x_max = ((x + r) <= w) ? (x + r): w;
      int y_max = ((y + r) <= h) ? (y + r): h;
      std::vector<float> norm_w(channels, 0.0);
      std::vector<float> val(channels, 0.0);

      // printf("Aqui\n");
      for (int j = y_min; j < y_max; j++){
      for (int i = x_min; i < x_max; i++){
        // printf("(%d,%d) (%d,%d)\n", x,y, i,j);
        float sp_w = exp(-(pow(i-x, 2) / (2 * pow(sigma_x, 2)) + pow(j-y, 2) / (2 * pow(sigma_y, 2)))); 
        for (int l = 0; l < channels; l ++){
          // printf("sp_w: %f, l:%d\n", sp_w, l);
          float weight = sp_w* exp(-pow(source(x,y,l)-source(i,j,l),2) / (2 * pow(sigma_c, 2)));
          norm_w[l] += weight;
          val[l] += source(i,j,l)* weight; 
        }
      }
      }
        
      for (int l = 0; l < channels; l++){
        val[l] /= norm_w[l];
        output(x,y,l) = val[l];
      } 

    }
    }
  }

  return output;
}