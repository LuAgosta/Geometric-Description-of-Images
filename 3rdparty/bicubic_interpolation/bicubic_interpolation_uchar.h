// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#ifndef BICUBIC_INTERPOLATION_UCHAR_H
#define BICUBIC_INTERPOLATION_UCHAR_H

#include <stdbool.h>

/**
  *
  * Compute the bicubic interpolation of a point in an image.
  * Detect if the point goes outside the image domain.
  *
**/
unsigned char bicubic_interpolation_uchar_at(
	const unsigned char *input, //image to be interpolated
	const unsigned char  uu,    //x component of the vector field
	const unsigned char  vv,    //y component of the vector field
	const int    nx,    //image width
	const int    ny,    //image height
  const int    pd,
	bool         border_out //if true, return zero outside the region
);


/**
  *
  * Compute the bicubic interpolation of an image.
  *
**/
void bicubic_interpolation_uchar_warp(
	const unsigned char *input,     // image to be warped
	const unsigned char *u,         // x component of the vector field
	const unsigned char *v,         // y component of the vector field
	unsigned char       *output,    // image warped with bicubic interpolation
	const int    nx,        // image width
	const int    ny,        // image height
  const int    pd,
  const int    l,
	bool         border_out // if true, put zeros outside the region
);

#endif//BICUBIC_INTERPOLATION_UCHAR_H
