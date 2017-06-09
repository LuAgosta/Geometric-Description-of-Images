// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#ifndef BICUBIC_INTERPOLATION_UCHAR_C
#define BICUBIC_INTERPOLATION_UCHAR_C

#include <stdbool.h>

#define BOUNDARY_CONDITION 0
//0 Neumann
//1 Periodic
//2 Symmetric

/**
  *
  * Neumann boundary condition test
  *
**/
static int neumann_bc(int x, int nx, bool *out)
{
	if(x < 0)
	{
	    x = 0;
	    *out = true;
	}
	else if (x >= nx)
	{
	    x = nx - 1;
	    *out = true;
	}

	return x;
}

/**
  *
  * Periodic boundary condition test
  *
**/
static int periodic_bc(int x, int nx, bool *out)
{
	if(x < 0)
	{
		const int n   = 1 - (int)(x/(nx+1));
		const int ixx = x + n * nx;

		x =   ixx% nx;
		*out = true;
	}
	else if(x >= nx)
	{
		x = x % nx;
		*out = true;
	}

	return x;
}


/**
  *
  * Symmetric boundary condition test
  *
**/
static int symmetric_bc(int x, int nx, bool *out)
{
	if(x < 0)
	{
		const int borde = nx - 1;
		const int xx = -x;
		const int n  = (int)(xx/borde) % 2;

		if ( n ) x = borde - ( xx % borde );
		else x = xx % borde;
		*out = true;
	}
	else if ( x >= nx )
	{
		const int borde = nx - 1;
		const int n = (int)(x/borde) % 2;

		if ( n ) x = borde - ( x % borde );
		else x = x % borde;
		*out = true;
	}

	return x;
}


/**
  *
  * Cubic interpolation in one dimension
  *
**/
static unsigned char cubic_interpolation_cell (
	unsigned char v[4],  //interpolation points
	unsigned char x      //point to be interpolated
)
{
	return  v[1] + 0.5 * x * (v[2] - v[0] +
		x * (2.0 *  v[0] - 5.0 * v[1] + 4.0 * v[2] - v[3] +
		x * (3.0 * (v[1] - v[2]) + v[3] - v[0])));
}


/**
  *
  * Bicubic interpolation in two dimensions
  *
**/
static unsigned char bicubic_interpolation_uchar_cell (
	unsigned char p[4][4], //array containing the interpolation points
	unsigned char x,       //x position to be interpolated
	unsigned char y        //y position to be interpolated
)
{
	unsigned char v[4];
	v[0] = cubic_interpolation_cell(p[0], y);
	v[1] = cubic_interpolation_cell(p[1], y);
	v[2] = cubic_interpolation_cell(p[2], y);
	v[3] = cubic_interpolation_cell(p[3], y);
	return cubic_interpolation_cell(v, x);
}

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
	const int 	 pd,   //Total Bytes number
	const int 	 l,    //Byte
	bool         border_out //if true, return zero outside the region
)
{
	const int sx = (uu < 0)? -1: 1;
	const int sy = (vv < 0)? -1: 1;

	int x, y, mx, my, dx, dy, ddx, ddy;
	bool out[1] = {false};

	//apply the corresponding boundary conditions
	switch(BOUNDARY_CONDITION) {

		case 0: x   = neumann_bc((int) uu, nx, out);
			y   = neumann_bc((int) vv, ny, out);
			mx  = neumann_bc((int) uu - sx, nx, out);
			my  = neumann_bc((int) vv - sx, ny, out);
			dx  = neumann_bc((int) uu + sx, nx, out);
			dy  = neumann_bc((int) vv + sy, ny, out);
			ddx = neumann_bc((int) uu + 2*sx, nx, out);
			ddy = neumann_bc((int) vv + 2*sy, ny, out);
			break;

		case 1: x   = periodic_bc((int) uu, nx, out);
			y   = periodic_bc((int) vv, ny, out);
			mx  = periodic_bc((int) uu - sx, nx, out);
			my  = periodic_bc((int) vv - sx, ny, out);
			dx  = periodic_bc((int) uu + sx, nx, out);
			dy  = periodic_bc((int) vv + sy, ny, out);
			ddx = periodic_bc((int) uu + 2*sx, nx, out);
			ddy = periodic_bc((int) vv + 2*sy, ny, out);
			break;

		case 2: x   = symmetric_bc((int) uu, nx, out);
			y   = symmetric_bc((int) vv, ny, out);
			mx  = symmetric_bc((int) uu - sx, nx, out);
			my  = symmetric_bc((int) vv - sx, ny, out);
			dx  = symmetric_bc((int) uu + sx, nx, out);
			dy  = symmetric_bc((int) vv + sy, ny, out);
			ddx = symmetric_bc((int) uu + 2*sx, nx, out);
			ddy = symmetric_bc((int) vv + 2*sy, ny, out);
			break;

		default:x   = neumann_bc((int) uu, nx, out);
			y   = neumann_bc((int) vv, ny, out);
			mx  = neumann_bc((int) uu - sx, nx, out);
			my  = neumann_bc((int) vv - sx, ny, out);
			dx  = neumann_bc((int) uu + sx, nx, out);
			dy  = neumann_bc((int) vv + sy, ny, out);
			ddx = neumann_bc((int) uu + 2*sx, nx, out);
			ddy = neumann_bc((int) vv + 2*sy, ny, out);
			break;
	}

	if(*out && border_out)
		return 0.0;

	else
	{
		//obtain the interpolation points of the image
		const unsigned char p11 = input[(mx  + nx * my)*pd + l];
		const unsigned char p12 = input[(x   + nx * my)*pd + l];
		const unsigned char p13 = input[(dx  + nx * my)*pd + l];
		const unsigned char p14 = input[(ddx + nx * my)*pd + l];

		const unsigned char p21 = input[(mx  + nx * y)*pd + l];
		const unsigned char p22 = input[(x   + nx * y)*pd + l];
		const unsigned char p23 = input[(dx  + nx * y)*pd + l];
		const unsigned char p24 = input[(ddx + nx * y)*pd + l];

		const unsigned char p31 = input[(mx  + nx * dy)*pd + l];
		const unsigned char p32 = input[(x   + nx * dy)*pd + l];
		const unsigned char p33 = input[(dx  + nx * dy)*pd + l];
		const unsigned char p34 = input[(ddx + nx * dy)*pd + l];

		const unsigned char p41 = input[(mx  + nx * ddy)*pd + l];
		const unsigned char p42 = input[(x   + nx * ddy)*pd + l];
		const unsigned char p43 = input[(dx  + nx * ddy)*pd + l];
		const unsigned char p44 = input[(ddx + nx * ddy)*pd + l];

		//create array
		unsigned char pol[4][4] = {
			{p11, p21, p31, p41},
			{p12, p22, p32, p42},
			{p13, p23, p33, p43},
			{p14, p24, p34, p44}
		};

		//return interpolation
		return bicubic_interpolation_uchar_cell(pol, uu-x, vv-y);
	}
}


/**
  *
  * Compute the bicubic interpolation of an image.
  *
**/
void bicubic_interpolation_uchar_warp(
	const unsigned char *input,     // image to be warped
	const float *u,         // x component of the vector field
	const float *v,         // y component of the vector field
	unsigned char       *output,    // image warped with bicubic interpolation
	const int    nx,        // image width
	const int    ny,        // image height
	const int    pd,				//Bytes number
	bool         border_out // if true, put zeros outside the region
)
{
	#pragma omp parallel for
	for(int i = 0; i < ny; i++)
		for(int j = 0; j < nx; j++)
			for (int l = 0; l < pd; l++){
				const int   p  = (i * nx + j)*pd + l;
				const float uu = (float) (j + u[p]);
				const float vv = (float) (i + v[p]);
				// obtain the bicubic interpolation at position (uu, vv)
				output[p] = bicubic_interpolation_uchar_at(input,
						uu, vv, nx, ny, pd, l,  border_out);
			}
}


#endif//BICUBIC_INTERPOLATION_UCHAR_C
