/**
 * Copyright (C) 2016, Roberto P.Palomares <roberto.palomares@upf.edu>
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD
 * License. You should have received a copy of this license along
 * this program. If not, see
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#ifndef FLST_UTILITY_H_
#define FLST_UTILITY_H_


extern "C" {
#include "flst.h"
}
#include <algorithm>
#include <vector>
#include <map>
#include <utility>
#include <set>

#include "image.h"
#include "gaussian_weights.h"
#include "bilateral_filter.h"
//TODO:Erase after debug
#include <stdio.h>

#include "io_utility.h"

/**
 * Contains methods to work with a tree shape of an image.
 * Acts as a proxy for Monase code
 */
class FLSTUtility
{
public:


  FLSTUtility(int pMinArea);
  FLSTUtility(int pMinArea, int radius);
  ~FLSTUtility();
  void calculate(FixedImage<float> source);
  FixedImage<float> get_cc_image_flow();

  //General getters and Setters
  int get_index_parent(int index);
  int get_index_first_child(int index);
  int get_index_next_sibling(int index);
  int get_index_sibling_list(int index);
  
  shape * get_shape(int index);
  shape * get_parent(int index);
  shape * get_child(int index);
  shape * get_next_sibling(int index);
  shape * get_sibling_list(int index);


  //Getters for the minimum bounding box
  int get_x_min_cc(int index, int cc);
  int get_y_min_cc(int index, int cc);
  int get_x_max_cc(int index, int cc);
  int get_y_max_cc(int index, int cc);


  //Getters and Setters for Shapes
  int get_pMin_Area();
  void set_pMin_Area(int pMinArea);
  int get_number_of_shapes();
  int get_region(int index);
  int get_area(int index);
  int get_tp_coord_x(int index, int pos);
  int get_tp_coord_y(int index, int pos);
  point_plane* get_tp_list(int index);

  //Getters and Setters for the connected component.
  int get_number_of_cc(int index);
  int get_number_elements_cc(int index,int cc);
  int get_tp_coord_x_cc(int index, int cc, int pos);
  int get_tp_coord_y_cc(int index, int cc, int pos);
  std::vector<float> get_of_cc(int index, int cc);
  void set_of_cc(std::vector<float> u, int index, int cc);
  std::vector<std::pair<float,float>> get_of_candidate_childs(int index);
  std::vector<std::pair<float,float>> get_of_candidate_parent(int index);
  std::vector<std::pair<float,float>> get_of_candidate_siblings(int index);
  std::vector<std::pair<float,float>> get_of_candidate_topological_neighbors(int index, int cc);
  std::vector<int> get_index_cc_of_pixel(int x, int y);

  //Check behaviour of flst.h (at the bottom of the .cpp file)
  Image<float> test_library();
  Image<float> test_library(FixedImage<float> source, int pMinArea);
  void draw_cc_shape(int index);
  void draw_topological_neighbors_cc(int index);
  void draw_all_cc();
  void test_region_pixels(int minArea);
  void test_at_least_one_pixel();
  void test_shape_without_holes();
  void test_no_repetition();


  //Total Variation, global and local, for an auxiliary image  /*Lucie*/
  float global_variation(Fimage auxImage) ;
  float local_variation(Fimage auxImage, int index, int cc) ;




private:
   int _pMinArea;
  shapes *_pTree;
  shape  **_root_to_leave;
  std::vector<std::pair<int,int>> _cc_pixel_map;
  Image<float> gaussian_weights_min_cc;
  int call_parent;
  int call_child ; 

  //Bounding box
  int _r_bb;
  float _sigma_c;

  //Create an order to move along the tree (pre-order) and set the index
  void topological_order();
  //Establish the number of private pixels of each region.
  void shape_without_holes();
  //Determine the cc for each shape.
  void create_conected_component_tree();
  //Determine the topological  neigbors (in the cc sense) for each cc.
  void create_topological_neigbors_list();
  //Create a correspondence map between pixels and their cc
  void create_cc_pixel_map();
  //Create mininum bounding box to contain each cc.
  void create_min_bounding_box_cc();
  //Create extended bounding box based on a radius (r)
  void create_extended_bounding_box_cc(int r);
  //Create gaussian weight (spatial and color)
  void create_gaussian_weigh_bounding_box(FixedImage<float> source,
                                                     float sigma_c);
  void create_gaussian_weigh_bounding_box_bilateral(FixedImage<float> source,
                                                  float sigma_c);
  //Liberate all memory reserved for the bounding box
  void liberate_bounding_box();

  //Produce the bilateral filter image
void create_bilateral_image(FixedImage<float> source);


  //Auxiliar functions to obtain cc. (Coloma code)
  void barOmega(int xO, int yO, int etiquetaCC,
                        double hole_graylev, int *original,
                        int *Omegah, int dx,int dy);
unsigned char hole(int xp, int yp, int etiquetaCC, 
      double ming, double maxg, int  *Image, int *OmegaH, unsigned char connecty8, 
      int dx, int dy);

Fimage create_auxImage_parent(int index , int cc)  ; /*Lucie*/
Fimage create_auxImage_child(int index, int cc)  ; /*Lucie*/
Fimage prune_cc() ; /*Lucie*/
//void prune_cc() ; /*Lucie*/
void prune_cc(Fimage im_out) ; /*Lucie*/
void modify_tree() ; /*Lucie*/
void create_extended_bounding_box_cc(int r, int index, int cc) ;  /*Lucie*/
Fimage create_patch_child(Fimage im_out, int index, int cc) ;
Fimage create_patch_parent(Fimage im_out, int index, int cc) ;
  //////////////////////////


  Fimage FixedImage_to_Fimage(FixedImage<float> source);
  FixedImage<float> Fimage_to_FixedImage(Fimage source);

  void getminmax(float *min,  float *max,  float *x,  int n);
  void image_normalization(float *a, float *an, int size);

    struct point_comparator
  {
    inline bool operator()(struct point_plane e1, struct point_plane e2)
    {
      if (e1.x != e2.x){
       return e1.x < e2.x;
      }
      return e1.y < e2.y;
    }
  };

  struct equal_comparator
  {
    inline bool operator()(struct point_plane e1, struct point_plane e2)
    {
      return ((e1.x == e2.x) && (e1.y == e2.y));
    }
  };


  struct shape_comparator
  {
    inline bool operator()(struct shape* e1, struct shape* e2)
    {
      if (e1->parent != e2->parent){
       return e1->parent < e2->parent;
      }
      return e1->child < e2->child;
    }
  };

  struct equal_shape_comparator
  {
    inline bool operator()(struct shape* e1, struct shape* e2)
    {
      return ((e1->parent == e2->parent) && (e1->child == e2->child)
        && (e1->next_sibling == e2->next_sibling));
    }
  };
};

#endif /* FLST_UTILITY_H_ */
