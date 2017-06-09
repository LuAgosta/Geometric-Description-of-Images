#include <stdio.h>
#include <stdlib.h>

#include "iio.h"
#include "flst.h"

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

  Fimage *imRGB = (Fimage*) calloc(3, sizeof(Fimage));
  imRGB[0] = mw_new_fimage();
  imRGB[1] = mw_new_fimage();
  imRGB[2] = mw_new_fimage();

  unsigned int i, ie;
  mw_change_fimage(imRGB[0], h, w);
  mw_change_fimage(imRGB[1], h, w);
  mw_change_fimage(imRGB[2], h, w);
  unsigned int nPixels = w * h;
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
    Fimage my_image = readImage(argv[1]);
    printf("\t* Image size: %dx%d\n", (int) my_image->ncol, (int) my_image->nrow);

    if (my_image) {

      Fimage im_out = mw_new_fimage();
      Shapes pTree = mw_new_shapes();
      int pMinArea = atoi(argv[3]);

      printf("\n* Apply FLST with parameter %d.\n", pMinArea);
      flst_bilinear(&pMinArea, my_image, pTree);
      // Alternatively
      // flst(&pMinArea, my_image, pTree);
      printf("* Reconstruction.\n");
      flst_reconstruct(pTree, im_out);
      printf("* Save the image as: %s.\n", argv[2]);
      saveImage(im_out, argv[2]);


      // FUNCTION TV : Does not work
      printf("\n* Tv function does not work.\n");
      // float pScale = 5, pQuantizationLevel = 1.0f;
      // int pQuantizationCurve = 4;
      // flstb_tv(&pScale, &pQuantizationLevel, &pQuantizationCurve, pTree, im_out);

      // FUNCTION QUANTIZE
      printf("\n* Test quantization function.\n");
      Fsignal pLevels = mw_new_fsignal();
      pLevels = mw_change_fsignal(pLevels, 5);
      pLevels->values[0] = 0;
      pLevels->values[1] = 64;
      pLevels->values[2] = 128;
      pLevels->values[3] = 192;
      pLevels->values[4] = 226;
      float pOffset = 12.5f, pStep = 50.0f;
      Shapes pQuantizedTree = mw_new_shapes();
      // Choose one of these two possibilities
      flstb_quantize(NULL, &pOffset, &pStep, pTree, pQuantizedTree);
      // flstb_quantize(pLevels, NULL, NULL, pTree, pQuantizedTree);
      printf("* Reconstruction.\n");
      flst_reconstruct(pQuantizedTree, im_out);
      printf("* Save the quantized image as: %s.\n", "imq.png");
      saveImage(im_out, "imq.png");

      Flist pBoundary = mw_new_flist();
      flst_boundary(pTree, pTree->the_shapes->child, pBoundary);


      printf("\n* Free memory\n");
      mw_delete_fimage(my_image);
      mw_delete_fimage(im_out);
      
    } else  {
      printf("Cannot read the images.\n");
      exit(1);
    }

  }
  return 0;
}
