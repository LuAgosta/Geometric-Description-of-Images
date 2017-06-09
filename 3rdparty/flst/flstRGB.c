#include <stdio.h>
#include <stdlib.h>

#include "iio.h"
#include "flst.h"
Fimage *mw_new_fimageRGB () {
  Fimage *imRGB = (Fimage*) calloc(3, sizeof(Fimage));
  imRGB[0] = mw_new_fimage();
  imRGB[1] = mw_new_fimage();
  imRGB[2] = mw_new_fimage();
  return imRGB;
}

void mw_change_fimageRGB (Fimage* imRGB, int h, int w) {
  imRGB[0] = mw_change_fimage(imRGB[0], h, w);
  imRGB[1] = mw_change_fimage(imRGB[1], h, w);
  imRGB[2] = mw_change_fimage(imRGB[2], h, w);
}

void mw_delete_fimageRGB (Fimage* imRGB) {
  mw_delete_fimage(imRGB[0]);
  mw_delete_fimage(imRGB[1]);
  mw_delete_fimage(imRGB[2]);
  free(imRGB);
}

Shapes* mw_new_shapesRGB () {
  Shapes* pTreeRGB = (Shapes*) calloc(3, sizeof(Shapes));
  pTreeRGB[0] = mw_new_shapes();
  pTreeRGB[1] = mw_new_shapes();
  pTreeRGB[2] = mw_new_shapes();
  return pTreeRGB;
}

void mw_delete_shapesRGB (Shapes* pTreeRGB) {
  mw_delete_shapes(pTreeRGB[0]);
  mw_delete_shapes(pTreeRGB[1]);
  mw_delete_shapes(pTreeRGB[2]);
  free(pTreeRGB);
}


void flst_bilinearRGB(int *pMinArea, Fimage *imRGB, Shapes *pTreeRGB) {
  flst_bilinear(pMinArea, imRGB[0], pTreeRGB[0]);
  flst_bilinear(pMinArea, imRGB[1], pTreeRGB[1]);
  flst_bilinear(pMinArea, imRGB[2], pTreeRGB[2]);
}
                                              

void flst_reconstructRGB(Shapes *pTree, Fimage *pFloatImageOutput) {
  flst_reconstruct(pTree[0], pFloatImageOutput[0]);
  flst_reconstruct(pTree[1], pFloatImageOutput[1]);
  flst_reconstruct(pTree[2], pFloatImageOutput[2]);
}







float *reshapeRGBImage(float *imin, unsigned int w, unsigned int h) {
  
  float *imout = (float *) calloc(3 * w * h, sizeof(float));
  unsigned int nPixels = w * h;

  float *Rin = imin;
  float *Gin = imin + 1;
  float *Bin = imin + 2;

  float *Rout = imout;
  float *Gout = imout + nPixels;
  float *Bout = imout + 2 * nPixels;

  float *ie = Rin + 3 * nPixels;
  for (; Rin < ie; Rin += 3, Gin += 3, Bin += 3, Rout++, Gout++, Bout++) {
    *Rout = *Rin;
    *Gout = *Gin;
    *Bout = *Bin;
  }
  free(imin);
  return imout;
}

float *unreshapeRGBImage(float *imin, unsigned int w, unsigned int h) {
  float *imout = (float *) calloc(3 * w * h, sizeof(float));
  unsigned int nPixels = w * h;
  float *Rin = imin, *Gin = imin + 1 * nPixels, *Bin = imin + 2 * nPixels;
  float *Rout = imout, *Gout = imout + 1, *Bout = imout + 2;
  float *ie = Rin + nPixels;
  for (; Rin < ie; Rin++, Gin++, Bin++, Rout += 3, Gout += 3, Bout += 3) {
    *Rout = *Rin;
    *Gout = *Gin;
    *Bout = *Bin;
  }
  free(imin);
  return imout;
}

Fimage readImage(const char* name) {
  int w, h, c;
  
  float *im = iio_read_image_float_vec(name, &w, &h, &c);
  if (c != 3) {
    printf("Function readImage: image must have 3 channels (RGB) not %d.\n", c);
    exit(1);
  }
  im = reshapeRGBImage(im, w, h);

  Fimage my_image = mw_new_fimage();
  unsigned int i, ie;
  mw_change_fimage(my_image, h, w);
  unsigned int nPixels = w * h;
  for (i = 0, ie = w * h; i < ie; i++) {
    my_image->gray[i] = (im[i] + im[i + nPixels] + im[i + 2 * nPixels]) / (3.0f);
  }
  return my_image;
}

Fimage* readImageRGB(const char* name) {
  int w, h, c;
  
  float *im = iio_read_image_float_vec(name, &w, &h, &c);
  if (c != 3) {
    printf("Function readImage: image must have 3 channels (RGB) not %d.\n", c);
    exit(1);
  }
  im = reshapeRGBImage(im, w, h);

  Fimage *imRGB = mw_new_fimageRGB();
  mw_change_fimageRGB(imRGB, h, w);

  unsigned int i, ie, nPixels = w * h;
  for (i = 0, ie = w * h; i < ie; i++) {
    imRGB[0]->gray[i] = im[i];
    imRGB[1]->gray[i] = im[i + nPixels];
    imRGB[2]->gray[i] = im[i + 2 * nPixels];
  }
  return imRGB;
}

void saveImage(Fimage my_image, const char *name) {
  unsigned int w = my_image->ncol, h = my_image->nrow;
  unsigned int nPixels = w * h;
  float *im = (float *) malloc(nPixels * 3 * sizeof(float));
  unsigned int i; 
  for (i = 0; i < nPixels; i++) {
    im[i] = (unsigned int) (my_image->gray[i]);  
    im[i + nPixels] = (unsigned int) (my_image->gray[i]);  
    im[i + 2 * nPixels] = (unsigned int) (my_image->gray[i]);  
  }
  im = unreshapeRGBImage(im, w, h);
  char fname[100];
  sprintf(fname, "%s", name);
  iio_save_image_float_vec(fname, im, w, h, 3);
  free(im);
}

void saveImageRGB (Fimage *imRGB, const char *name) {
  unsigned int w = imRGB[0]->ncol, h = imRGB[0]->nrow;
  unsigned int nPixels = w * h;
  float *im = (float *) malloc(nPixels * 3 * sizeof(float));
  unsigned int i; 
  for (i = 0; i < nPixels; i++) {
    im[i] = (unsigned int) (imRGB[0]->gray[i]);  
    im[i + nPixels] = (unsigned int) (imRGB[1]->gray[i]);  
    im[i + 2 * nPixels] = (unsigned int) (imRGB[2]->gray[i]);  
  }
  im = unreshapeRGBImage(im, w, h);
  char fname[100];
  sprintf(fname, "%s", name);
  iio_save_image_float_vec(fname, im, w, h, 3);
  free(im);
}

int main(int argc, char *argv[]) {

  if (argc < 4) {

    printf("* Usage : flst <input> <output> <min_area>\n");

  } else  {
    printf("* Read the image: %s\n", argv[1]);
    Fimage *imRGB = readImageRGB(argv[1]);
    printf("\t* Image size: %dx%d\n", (int) imRGB[0]->ncol, (int) imRGB[0]->nrow);

    if (imRGB) {

      Fimage *imRGBOut = mw_new_fimageRGB();
      Shapes *pTreeRGB = mw_new_shapesRGB();
      int pMinArea = atoi(argv[3]);

      printf("\n* Apply FLST with parameter %d.\n", pMinArea);
      flst_bilinearRGB(&pMinArea, imRGB, pTreeRGB);


      // FUNCTION QUANTIZE
      printf("\n* Test quantization function.\n");
      Fsignal pLevels = mw_new_fsignal();
      pLevels = mw_change_fsignal(pLevels, 5);
      pLevels->values[0] = 0;
      pLevels->values[1] = 64;
      pLevels->values[2] = 128;
      pLevels->values[3] = 192;
      pLevels->values[4] = 226;
      float pOffset = 5.0f, pStep = 10.0f;
      Shapes *pQuantizedTree = mw_new_shapesRGB();
      // Choose one of these two possibilities
      flstb_quantize(NULL, &pOffset, &pStep, pTreeRGB[0], pQuantizedTree[0]);
      flstb_quantize(NULL, &pOffset, &pStep, pTreeRGB[1], pQuantizedTree[1]);
      flstb_quantize(NULL, &pOffset, &pStep, pTreeRGB[2], pQuantizedTree[2]);

      printf("* Reconstruction.\n");
      flst_reconstructRGB(pQuantizedTree, imRGBOut);
      printf("* Save the quantized image as: %s.\n", "imq.png");
      saveImageRGB(imRGBOut, "imq.png");



      printf("* Reconstruction.\n");
      flst_reconstructRGB(pTreeRGB, imRGBOut);

      printf("* Save the image as: %s.\n", argv[2]);
      saveImageRGB(imRGBOut, argv[2]);

      printf("\n* Free memory\n");
      mw_delete_shapesRGB(pTreeRGB);
      mw_delete_fimageRGB(imRGB);
      mw_delete_fimageRGB(imRGBOut);
     
    } else  {
      printf("Cannot read the images.\n");
      exit(1);
    }

  }
  return 0;
}
