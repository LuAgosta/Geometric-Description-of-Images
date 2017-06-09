/**
 * Copyright (C) 2016, Roberto P.Palomares <roberto.palomares@upf.edu>
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD
 * License. You should have received a copy of this license along
 * this program. If not, see
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include <string>
#include <ctime>

// #include "image_inpainting.h"
#include "io_utility.h"
#include "flst_utility.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdexcept>

// NOTE: define DBG_OUTPUT in image_inpainting.h header to turn on the output after each iteration and other debug images output

using namespace std;

/**
 * @param c pointer to original argc
 * @param v pointer to original argv
 * @param name option name after hyphen
 * @param default_value default value (if NULL, the option takes no argument)
 */
static const char *pick_option(int *c, char ***v, const char *name, const char *default_value)
{
  int argc = *c;
  char **argv = *v;
  int id = default_value ? 1 : 0;
  for (int i = 0; i < argc - id; i++) {
    if (argv[i][0] == '-' && 0 == strcmp(argv[i] + 1, name)) {
      char *r = argv[i + id] + 1 - id;
      *c -= (id + 1);
      for (int j = i; j < argc - id; j++) {
        (*v)[j] = (*v)[j + id + 1];
      }

      return r;
    }
  }

  return default_value;
}


void  add_margin(Mask first_mask, Mask second_mask, int margin)
{
  for (int i = 0; i < margin; i++) {
    for (uint x = 0; x < first_mask.get_size_x(); x++) {
      first_mask.unmask(x, i);
      first_mask.unmask(x, first_mask.get_size_y() - i - 1);
      second_mask.unmask(x, i);
      second_mask.unmask(x, second_mask.get_size_y() - i - 1);
    }
    for (uint y = 0; y < first_mask.get_size_y(); y++) {
      first_mask.unmask(i, y);
      first_mask.unmask(first_mask.get_size_x() - i - 1, y);
      second_mask.unmask(i, y);
      second_mask.unmask(second_mask.get_size_x() - i - 1, y);
    }
  }
}



int main(int argc, char *argv[])
{


  // get parameters from command line
  int pMinArea           = atoi(pick_option(&argc, &argv, "area"  , "20"));      
  

  if (argc < 3) {
    // display usage message and quit
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s a output [OPTIONS]\n\n", argv[0]);
    fprintf(stderr, "Available options are:\n");
    fprintf(stderr, " -area  \tpatch side (%d)\n", pMinArea);
    return 1;
  }

  // get image paths from command line (or use default values)
  int i=1;
  string a_name  = (argc>i) ? argv[i]: "a.png"; i++;
  string output_name = (argc>i) ? argv[i]: "out.png"; i++;
  //printf("%s %s %s\n", input_name.c_str(), mask_name.c_str(), output_name.c_str());


  //Creation of Images channel = 3 /*Lucie*/
  // FixedImage<float> *newImage = new FixedImage<float>(500, 500, (uint)3); 
  // newImage->draw_rect (0 , 0, 499, 499,0, 50) ; 
  // newImage->draw_rect (100 , 100, 300, 300,1, 200) ; 
  // newImage->draw_rect (150, 90, 100, 9,2,  100) ; 
  // newImage->draw_rect (401, 150, 10, 50,0, 100) ; 
  // IOUtility::write_rgb_image(a_name, *newImage);
  // printf("number of channels of newImage: %d \n" , newImage->get_number_of_channels()) ;

  //Creation of Images channel = 1 /*Lucie*/
  // FixedImage<float> *newImage = new FixedImage<float>(500, 500,(float)50); 
  // newImage->draw_rect (100 , 100, 300, 300, 200) ; 
  // newImage->draw_rect (150, 90, 100, 9, 100) ; 
  // newImage->draw_rect (401, 150, 10, 50, 100) ; 
  // IOUtility::write_mono_image(a_name, *newImage);
  // printf("number of channels of newImage: %d \n" , newImage->get_number_of_channels()) ;


#ifdef DBG_OUTPUT
  // define prefix for debug output
  IOUtility::set_prefix("dbg/");
#endif

  // read image and mask
  Image<float> input_a = IOUtility::rgb_to_lab(IOUtility::read_rgb_image(a_name));
 
  // display the stored image a /*Lucie*/
  printf("number of channels of image converted in Fimage : %d \n" , input_a.get_number_of_channels()) ; /*Lucie*/
  /*channel = 3 always ! */
  IOUtility::write_rgb_image("Image used " + a_name , input_a) ; /*Lucie*/
  

  //Prueba funcionamiento FLST
  FLSTUtility *flst_utility = new FLSTUtility(pMinArea);
  flst_utility->calculate(input_a);
  Image<float> prueba_a = flst_utility->test_library();
  IOUtility::write_mono_image(output_name, prueba_a);
  printf("number of channels of output : %d \n" , prueba_a.get_number_of_channels()) ; /*Lucie*/


  delete flst_utility;
  //delete newImage ; 
}
