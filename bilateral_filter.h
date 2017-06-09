/**
 * Copyright (C) 2016, Roberto P.Palomares <roberto.palomares@upf.edu>
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD
 * License. You should have received a copy of this license along
 * this program. If not, see
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#ifndef BILATERAL_FILTER_H_
#define BILATERAL_FILTER_H_

#include <vector>
#include "image.h"

/**
 * Contains methods for calculating discrete Gaussian kernel.
 */
namespace BilateralFilter {
  Image<float> calculate(FixedImage<float> source,
                         int r,
                         int sigma_x,
                         int sigma_y,
                         int sigma_c,
                         int iter);
} // namespace BilateralFilter


#endif /* BILATERAL_FILTER_H_ */