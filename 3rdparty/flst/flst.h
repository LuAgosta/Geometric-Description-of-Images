#include "mw_error.h"
#include "shape.h"
#include "fimage.h"
#include "fsignal.h"
#include "list.h"


void flst(int* pMinArea, Fimage pImageInput, Shapes pTree);
void flst_bilinear(int* pMinArea, Fimage pImageInput, Shapes pTree);

void flst_reconstruct(Shapes pTree, Fimage pFloatImageOutput);

void flst_pixels(Shapes pTree);
Flist flst_boundary(Shapes pTree, Shape pShape, Flist pBoundary);

void flstb_dual(Shapes pTree, Shapes pDualTree);
void flstb_dualchain(Shapes, Shape, Flist, char*);
Flist flstb_boundary(int*, Fimage, Shapes, Shape, Flist, Flist, char*);

void flstb_quantize(Fsignal, float *, float *, Shapes, Shapes);



