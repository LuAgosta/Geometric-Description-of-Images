/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   shapes.h
   
   Vers. 1.2
   (C) 1999-2001 Pascal Monasse, Frederic Guichard, Jacques Froment.

   Internal Input/Output for the following FLST (Fast Level Set Transform)
   structures
     Point_plane
     Shape
     Shapes

Changes :
0.0: First version
0.1: minor modifications = some variable names (L.Moisan)
1.0: include "string_size.h" removed.
1.1: added data_size field in Shape and Shapes (L.Moisan), interpolation added.
1.2: no more modes for flst (L.Moisan)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~  This file is part of the MegaWave2 system library ~~~~~~~~~~~~~~~
MegaWave2 is a "soft-publication" for the scientific community. It has
been developed for research purposes and it comes without any warranty.
The last version is available at http://www.cmla.ens-cachan.fr/Cmla/Megawave
CMLA, Ecole Normale Superieure de Cachan, 61 av. du President Wilson,
      94235 Cachan cedex, France. Email: megawave@cmla.ens-cachan.fr 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifndef shape_flg
#define shape_flg

#ifdef SunOS
#include <sys/types.h>
#endif
#ifndef list_flg
#include "list.h"
#endif

#define mw_cmtsize  255 /* Size of the comments field */
#define mw_namesize 255 /* Size of the comments field */

/* Point in the plane */

typedef struct point_plane {
  int x,y; /* Coordinates of the point */
} *Point_plane;

/////ROBERTO/////////
typedef struct of_value {
  float u,v; /* Coordinates of the point */
} *Of_value;

typedef struct tp_neighbor{
  int index,cc; /* Coordinates of the point */
} *Tp_neighbor;
/////////////

/* A shape  : a connected component of a level set, with filled holes */
typedef struct shape
{
  char inferior_type; /* Indicates if it is extracted from a superior 
			 or inferior level set */
  float value; /* Limiting gray-level of the level set */
  char open; /* Indicates if the shape meets the border of the image */
  int area; /* Area of the shape = area of the cc of level set 
	                         + areas of the holes */
  char removed; /* Indicates whether the shape exists or not */

  Point_plane pixels; /* The array of pixels contained in the shape */

  //////////////////////////////////////////////////////////////////////////////
  //ROBERTO 
  //Connected component.
  int region; /* Number of private pixels*/
  int index;  /*Index to acces directly to the shape*/
  int cc_n; /*Number of connected component*/
  Point_plane *cc_ini; /*Pixels contained in each cc*/
  int *cc_nelem; /*Array with the number of pixels of each cc.*/
  Of_value cc_of; /*Array with the optical flow value for each cc*/
  //ROBERTO - Bounding box 
  int *x_min; //Contain the left- top corner for each cc
  int *y_min; //Contain the left- top corner for each cc
  int *x_max; //Contain the right - bottom corner for each cc
  int *y_max; //Contain the right - bottom corner for each cc
  float **gaussian_weight;
  //Topological neighbors
  int *tp_n; /*Array with the topological neighbors for each cc*/
  Tp_neighbor* tp_list; //Array of arrays (index, cc) for each cc*/
  //////////////////////////////////////////////////////////////////////////////


  Flist boundary; /* The boundary curve defining the shape */

  /* Data to include it in a tree. It has a parent (the smallest containing 
     shape), children (the largest contained shapes, whose first is pChild 
     and the others are its siblings), and siblings (the other children of 
     its parent) */
  struct shape *parent, *next_sibling, *child;

  int data_size;     /* size of data[] in bytes */
  void* data;        /* User defined field (saved). A pointer to something */

} *Shape;


/* A set of shapes (complete representation of an image) */
typedef struct shapes
{
  char cmt[mw_cmtsize];   /* Comments */
  char name[mw_namesize]; /* Name of the set */
  int nrow;               /* Number of rows (dy) of the image */
  int ncol;               /* Number of columns (dx) of the image */
  int interpolation;      /* Interpolation used for the level lines:
			     0=nearest neighbor, 1=bilinear */
  Shape the_shapes; /* Array of the shapes. 
		       The root of the tree is at index 0 */

  int nb_shapes; /* The number of shapes (the size of the array the_shapes) */

  /* Link between pixels and shapes */
  Shape *smallest_shape; /* An image giving for each pixel 
			    the smallest shape containing it */

  int data_size;     /* size of data[] in bytes */
  void* data;        /* User defined field (saved). A pointer to something */

} *Shapes;


/* Functions definition */

#ifdef __STDC__

Point_plane mw_new_point_plane(void);
Point_plane mw_change_point_plane(Point_plane);
Shape mw_new_shape(void);
Shape mw_change_shape(Shape);
void mw_delete_shape(Shape);
Shape mw_get_not_removed_shape(Shape);
Shape mw_get_parent_shape(Shape);
Shape mw_get_first_child_shape(Shape);
Shape mw_get_next_sibling_shape(Shape);
Shape mw_get_smallest_shape(Shapes, int, int);
Shapes mw_new_shapes(void);
Shapes mw_alloc_shapes(Shapes,int,int,float);
Shapes mw_change_shapes(Shapes,int,int,float);     
void mw_delete_shapes(Shapes);

#else

Point_plane mw_new_point_plane();
Point_plane mw_change_point_plane();
Shape mw_new_shape();
Shape mw_change_shape();
void mw_delete_shape();
Shape mw_get_not_removed_shape();
Shape mw_get_parent_shape();
Shape mw_get_first_child_shape();
Shape mw_get_next_sibling_shape();
Shape mw_get_smallest_shape();
Shapes mw_new_shapes();
Shapes mw_alloc_shapes();
Shapes mw_change_shapes();     
void mw_delete_shapes();

#endif

#endif

