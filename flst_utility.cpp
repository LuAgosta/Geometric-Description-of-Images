/**
 * Copyright (C) 2016, Roberto P.Palomares <roberto.palomares@upf.edu>
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD
 * License. You should have received a copy of this license along
 * this program. If not, see
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */
#include "flst_utility.h"


FLSTUtility::FLSTUtility(int pMinArea)
{
  _pMinArea = pMinArea;
  _pTree = mw_new_shapes();
  //TODO: Fixed value for the expanding bounding box.
  _r_bb = 5; 
  _sigma_c = 2;
}

FLSTUtility::FLSTUtility(int pMinArea, int radius)
{
  _pMinArea = pMinArea;
  if (radius >=0){
    _r_bb = radius;
  } else{
     _r_bb = 0;
  }
  //TODO:Sigma color
  _sigma_c = 2;
  _pTree = mw_new_shapes();
}

FLSTUtility::~FLSTUtility()
{
  if (_pTree != NULL){
    // Free Conected component structure
    for (int index = 0; index < get_number_of_shapes(); index++){
      shape *sh = get_shape(index);
      if (sh->cc_ini!=NULL)
        free(sh->cc_ini);
      if (sh->cc_nelem!=NULL)
        free(sh->cc_nelem);
    }
    // Free gaus weight.
    //Coloma: Comentado lo de la bounding box. No hace falta
    // liberate_bounding_box();
    mw_delete_shapes(_pTree);
  }
  free(_root_to_leave);
}



void FLSTUtility::calculate(FixedImage<float> source)
{
  if (_pMinArea < 0)
    _pMinArea = 1;
  Fimage input = FixedImage_to_Fimage(source);
  // give the global or local total variation of an Fimage                                         
  float TV = global_variation(input) ;
  printf("TV = %f \n", TV );
  flst(&_pMinArea, input, _pTree);
  

  flst_pixels(_pTree);
  topological_order();
  shape_without_holes();

  //Create conected component and compute their topological neighbors.
  create_conected_component_tree();
  create_cc_pixel_map();
  create_topological_neigbors_list();

  //Coloma: Esto creaba las bounding box. LO comento porque note hace falta
  //Bounding box function. Create a square patch around each cc.
  create_min_bounding_box_cc();
  // create_extended_bounding_box_cc(_r_bb);
  // create_gaussian_weigh_bounding_box(source, _sigma_c);
 
  //Chose the gray value of the pruned cc /*Lucie*/ 


  // Fimage auxImage_parent = create_auxImage_parent(3, 0) ; 
  // float TV_parent = global_variation(auxImage_parent) ;
  // printf("TV_parent = %f \n", TV_parent );
  // FixedImage<float> outputParent = Fimage_to_FixedImage(auxImage_parent);
  // IOUtility::write_mono_image("auxImageParent.png ", outputParent) ;


  // Fimage auxImage_child = create_auxImage_child(3 , 0) ;
  // float TV_child = global_variation(auxImage_child) ;
  // printf("TV_child = %f \n", TV_child ); 
  // FixedImage<float> outputChild = Fimage_to_FixedImage(auxImage_child);
  // IOUtility::write_mono_image("auxImageChild.png ", outputChild) ;


  // Fimage auxImage = prune_cc() ;   // diverge 
  // FixedImage<float> output = Fimage_to_FixedImage(auxImage);
  // IOUtility::write_mono_image("auxImage.png ", output) ;


  draw_all_cc();
  // create_bilateral_image(source);
  mw_delete_fimage(input);
  //mw_delete_fimage(auxImage);
  //mw_delete_fimage(auxImage_child);
  //mw_delete_fimage(auxImage_parent);
}


//Reconstruct the optical flow given by the tree shapes of an image.
FixedImage<float> FLSTUtility::get_cc_image_flow(){
  int h = _pTree->nrow;
  int w = _pTree->ncol;
  Image<float> flow(w, h, (uint)2);

  for (int index = 0; index < _pTree->nb_shapes; index++){
    int cc_n = get_number_of_cc(index);
    for (int l = 0; l < cc_n; l++){
      int         cc_nelem = get_number_elements_cc(index, l);
      std::vector<float> u = get_of_cc(index, l);
      for (int j = 0; j < cc_nelem; j++){
        const int x = get_tp_coord_x_cc(index, l, j);
        const int y = get_tp_coord_y_cc(index, l, j);
        flow(x, y, 0) = u[0];
        flow(x, y, 1) = u[1];
      }
    }
  }

  return flow;
}


//General getters and setters
/*General getters and setters*/
int FLSTUtility::get_index_parent(int index){
  if ((index <0) || (index >= _pTree->nb_shapes))
    return -1;
  if (_root_to_leave[index]->parent==NULL)
    return -1;
  return _root_to_leave[index]->parent->index;
}


int FLSTUtility::get_index_first_child(int index){
  if ((index <0) || (index >= _pTree->nb_shapes)
      || (_root_to_leave[index]->child==NULL))
  return -1;
  if (_root_to_leave[index]->child==NULL)
    return -1;
  return _root_to_leave[index]->child->index;
}

int FLSTUtility::get_index_next_sibling(int index){
  if ((index <0) || (index >= _pTree->nb_shapes))
    return -1;

  if (_root_to_leave[index]->next_sibling==NULL)
    return -1;

  return _root_to_leave[index]->next_sibling->index;
}

int FLSTUtility::get_index_sibling_list(int index){
  if ((index <0) || (index >= _pTree->nb_shapes))
  return -1;
  if (_root_to_leave[index]->parent==NULL)
    return -1;
  return _root_to_leave[index]->parent->child->index;
}


shape *FLSTUtility::get_shape(int index){
  if ((index <0) || (index >= _pTree->nb_shapes))
    return NULL;
  return _root_to_leave[index];
}

shape *FLSTUtility::get_parent(int index){
  if ((index <0) || (index >= _pTree->nb_shapes))
    return NULL;
  return _root_to_leave[index]->parent;
}

shape *FLSTUtility::get_child(int index){
  if ((index <0) || (index >= _pTree->nb_shapes))
    return NULL;
  return _root_to_leave[index]->child;
}

shape *FLSTUtility::get_next_sibling(int index){
  if ((index <0) || (index >= _pTree->nb_shapes))
    return NULL;
  return _root_to_leave[index]->next_sibling;
}

shape *FLSTUtility::get_sibling_list(int index){
  if ((index <0) || (index >= _pTree->nb_shapes))
    return NULL;

  if (get_parent(index)==NULL){//the root
    return get_child(index);
  }
  return _root_to_leave[index]->parent->child;
}

int FLSTUtility::get_x_min_cc(int index, int cc){
  return _root_to_leave[index]->x_min[cc];
}
int FLSTUtility::get_y_min_cc(int index, int cc){
  return _root_to_leave[index]->y_min[cc];
}
int FLSTUtility::get_x_max_cc(int index, int cc){
  return _root_to_leave[index]->x_max[cc];
}
int FLSTUtility::get_y_max_cc(int index, int cc){
  return _root_to_leave[index]->y_max[cc];
}


/*Setters and getters*/
int FLSTUtility::get_pMin_Area(){
  if (_pTree!=NULL)
    return _pMinArea;
  return 0;
}

void FLSTUtility::set_pMin_Area(int pMinArea){
  if (pMinArea> 0)
    _pMinArea = pMinArea;
}

int FLSTUtility::get_number_of_shapes(){
  if (_pTree!=NULL)
    return _pTree->nb_shapes;
  return 0;
}

int FLSTUtility::get_area(int index){
  if ((index <0) || (index > _pTree->nb_shapes))
    return -1;
  return _root_to_leave[index]->area;
}

int FLSTUtility::get_region(int index){
  if ((index <0) || (index > _pTree->nb_shapes))
    return -1;
  return _root_to_leave[index]->region;
}

int FLSTUtility::get_tp_coord_x(int index, int pos){
    if ((index <0) || (index >= _pTree->nb_shapes)
      || (pos <0 ) || (pos >= _root_to_leave[index]->region))
      return -1;
    return _root_to_leave[index]->pixels[pos].x;
}

int FLSTUtility::get_tp_coord_y(int index, int pos){
    if ((index <0) || (index >= _pTree->nb_shapes)
      || (pos <0 ) || (pos >= _root_to_leave[index]->region))
      return -1;
    return _root_to_leave[index]->pixels[pos].y;
}

point_plane* FLSTUtility::get_tp_list(int index){
  if ((index <0) || (index >= _pTree->nb_shapes))
    return NULL;
  return _root_to_leave[index]->pixels;
}

/*CC - Setters and getters*/
int FLSTUtility::get_number_of_cc(int index){
  if ((index <0) || (index >= _pTree->nb_shapes))
    return -1;
  return _root_to_leave[index]->cc_n;
}


int FLSTUtility::get_number_elements_cc(int index, int cc){
  if (index < get_number_of_shapes() && cc < get_number_of_cc(index)){
    return _root_to_leave[index]->cc_nelem[cc];
  }
  return -1;
}

int FLSTUtility::get_tp_coord_x_cc(int index, int cc, int pos){
    if ((index <0) || (index >= _pTree->nb_shapes)
      || (cc < 0 ) || (cc >= _root_to_leave[index]->cc_n)
      || (pos< 0 ) || (pos >= _root_to_leave[index]->cc_nelem[cc]))
      return -1;
    return _root_to_leave[index]->cc_ini[cc][pos].x;
}

int FLSTUtility::get_tp_coord_y_cc(int index, int cc, int pos){
    if ((index <0) || (index >= _pTree->nb_shapes)
      || (cc < 0 ) || (cc >= _root_to_leave[index]->cc_n)
      || (pos< 0 ) || (pos >= _root_to_leave[index]->cc_nelem[cc]))
      return -1;
    return _root_to_leave[index]->cc_ini[cc][pos].y;
}

std::vector<float> FLSTUtility::get_of_cc(int index, int cc){
  std::vector<float> u(2);
  if (index < get_number_of_shapes() && index >=0
      && cc < get_number_of_cc(index) && cc >=0){
    u[0] = _root_to_leave[index]->cc_of[cc].u;
    u[1] = _root_to_leave[index]->cc_of[cc].v;
    return u;
  } 
  return {};
}

void FLSTUtility::set_of_cc(std::vector<float> u, int index, int cc){
  if (index < get_number_of_shapes() && cc < get_number_of_cc(index)
      && index >=0 && cc >=0 && !u.empty()){
    _root_to_leave[index]->cc_of[cc].u = u[0];
    _root_to_leave[index]->cc_of[cc].v = u[1];
  } 
}

//Give a list with optical flow values from the childs
std::vector<std::pair<float,float>> FLSTUtility::get_of_candidate_childs(int index){
    std::vector< std::pair<float,float> > list;

    for (shape *sh = get_child(index); sh!= NULL; sh = sh->next_sibling){
      int cc_n = get_number_of_cc(sh->index);
      for (int l = 0; l < cc_n; l++){
        std::vector<float> u = get_of_cc(sh->index, l);
        std::pair<float,float> of(u[0], u[1]);
        list.push_back(of);
      }
    }
    return list;
}

//Give a list with all the optical flow from the parent
std::vector<std::pair<float,float>> FLSTUtility::get_of_candidate_parent(int index){
    std::vector< std::pair<float,float> > list;

    shape *sh = get_parent(index);
    if (sh!=NULL){
      int cc_n = get_number_of_cc(sh->index);
      for (int l = 0; l < cc_n; l++){
        std::vector<float> u = get_of_cc(sh->index, l);
        std::pair<float,float> of(u[0], u[1]);
        list.push_back(of);
       }
    }
    return list;
}

//Give a list with all the optical flow from the siblings
std::vector<std::pair<float,float>> FLSTUtility::get_of_candidate_siblings(int index){
    std::vector< std::pair<float,float> > list;

    shape *sh = get_sibling_list(index);
    shape *actual = get_shape(index);
    for (; sh!= NULL; sh = sh->next_sibling){
      if (sh!=actual){
        int cc_n = get_number_of_cc(sh->index);
        for (int l = 0; l < cc_n; l++){
          std::vector<float> u = get_of_cc(sh->index, l);
          std::pair<float,float> of(u[0], u[1]);
          list.push_back(of);
        }
      }
    }
    return list;
}

//Give a list with all the optical flow from the topological neighbors
std::vector<std::pair<float,float>> FLSTUtility::get_of_candidate_topological_neighbors(int index, int cc){
    std::vector< std::pair<float,float> > list;

    shape *sh = get_shape(index);
    int neig_n = sh->tp_n[cc];
    for (int i = 0; i < neig_n; i++){
      const int neig_index = sh->tp_list[cc][i].index;
      const int neig_cc = sh->tp_list[cc][i].cc;
      std::vector<float> u = get_of_cc(neig_index, neig_cc);
      std::pair<float,float> of(u[0], u[1]);
      list.push_back(of);
    }
    return list;
}


//Return a vector with (index, cc) from an specific pixel. It is obtained by linear search. It is nof efficient.
std::vector<int> FLSTUtility::get_index_cc_of_pixel(int x, int y){
  int h = _pTree->nrow;
  int w = _pTree->ncol;
  shape **smallest_shape = _pTree->smallest_shape;
  shape *sh;
  Image<float> flow(w, h, (uint)2);
  sh = smallest_shape[y*w + x];
  printf("Encuentra la shape\n");
  std::vector<int> id(2,-1);
  int index = sh->index;
  id[0] = index;
  int flag = 0;
  const int cc_n = get_number_of_cc(index);
  for (int l = 0; l < cc_n; l++){
    int cc_nelem = get_number_elements_cc(index,l);
    for (int j = 0; j < cc_nelem; j++){
      const int pos_x = get_tp_coord_x_cc(index, l, j);
      const int pos_y = get_tp_coord_y_cc(index, l, j);
      if ((pos_x == x) && (pos_y == y)){
        id[1] = l;
        flag++;
      }
    }
  }
  if (flag != 1){
    printf("Ha encontrado más de uno\n");
  }
  if ((id[0] <0) && (id[1] < 0)){
    printf("No lo ha encontrado :(\n");
  }
  return id;
}


/*Private*/
void FLSTUtility::topological_order(){
  if (!_pTree)
      return;

  // Shape root = _pTree->the_shapes;
  int size = _pTree->nb_shapes;

  _root_to_leave = (shape **) malloc(size*sizeof(shape *));
  for (int i= 0; i < size; i++){
    _root_to_leave[i] = &_pTree->the_shapes[i]; 
    (&_pTree->the_shapes[i])->index = i;
  }

}

//Establish the number of private pixels of each region.
void FLSTUtility::shape_without_holes(){
  shape *sh;
    int counter_01 = 0;
    int counter_20 = 0;
    for (int k = 0; k < _pTree->nb_shapes; k++){
      sh = &(_pTree->the_shapes[k]);
      point_plane *child_pixels;
      if (sh->child!=NULL){
        child_pixels = sh->child->pixels;
        point_plane *aux;
        int j = 0;
        aux = sh->pixels;
        while (aux!=child_pixels){
          aux++;
          j++;
        }
        sh->region = j;
      }else{
        sh->region = sh->area;
      }
      if (sh->region == 1){
        counter_01++;
      }
      if (sh->region >= 20){
        counter_20++;
      }
    }
    std::printf("Number of shapes: %d\n ",_pTree->nb_shapes);
    std::printf("Porcentaje (==1): %f\n", (counter_01*1.0)/_pTree->nb_shapes);
    std::printf("Porcentaje (>=20): %f\n", (counter_20*1.0)/_pTree->nb_shapes);

}


//Determine the cc of the tree.
void FLSTUtility::create_conected_component_tree(){

  int nb_shapes = get_number_of_shapes();
  int w = _pTree->ncol;
  int h = _pTree->nrow;

  //An auxiliary image initialize to 0.
  std::vector<int> labels_init(w*h,0);
  std::vector<int> labels_final(w*h,0);
  //Create primary label for each connected component inside of a region.
  for (int index = 0; index < nb_shapes; index++){
    int or_label = index +1;
    int nb_region = get_region(index);
    for (int l = 0; l < nb_region; l++){
      const int x =   get_tp_coord_x(index, l);
      const int y =   get_tp_coord_y(index, l);
      labels_init[y*w + x] = or_label;
    }
  }
  //TODO: Borrar printf
  printf("Imagen inicial con una etiqueta por region\n");

  //Extend the initial labels. Each region is divided by its connected component.
  int new_label = nb_shapes + 1; //
  std::vector<std::vector<int>> key_label(nb_shapes);
  for (int index = 0; index < nb_shapes; index++){
    int or_label = index +1;
    int nb_region = get_region(index);
    for (int l = 0; l < nb_region; l++){
      const int x =   get_tp_coord_x(index, l);
      const int y =   get_tp_coord_y(index, l);
      if (l == 0){ //First element to the region
        key_label[index].push_back(or_label);
        barOmega(x, y, or_label, or_label, 
                  &labels_init.front(), &labels_final.front(), w, h);
      }else if (labels_final[y*w + x]==0){ //If it is not zero, then it belongs to a cc.
        key_label[index].push_back(new_label);
        barOmega(x, y, new_label, or_label, 
                &labels_init.front(), &labels_final.front(), w, h);
        
        new_label ++;
      }    
    }
  }
  //TODO: Borrar printf
  printf("Imagen de etiquetas extendida\n");
  int total_CC = 0;
  int total_CC_mayores = 0;
  //Obtain the cc for each region
  for (int index = 0; index < nb_shapes; index++){
    int nb_region = get_region(index);
    std::map<int, std::vector<std::pair<int,int>>> map_label;

    // Create empty arrays.
    for(std::vector<int>::size_type i = 0; i != key_label[index].size(); i++) {
      std::vector<std::pair<int,int>> elem;
      map_label[key_label[index][i]] = elem;
    }

    //Insert all the elements of each connected compoment.
    for (int l = 0; l < nb_region; l++){
      const int x = get_tp_coord_x(index, l);
      const int y = get_tp_coord_y(index, l);
      map_label[labels_final[y*w + x]].push_back(std::make_pair(x,y));
    } 

    //Obtain the point plane of cc sortered
    std::vector<struct point_plane> pixels_cc;
    std::vector<int> cc_nelem;
    for(std::vector<int>::size_type i = 0; i != key_label[index].size(); i++) {
      std::vector<std::pair<int,int>> elem = map_label[key_label[index][i]];
      for (std::vector<std::pair<int,int>>::iterator it = elem.begin(); it!= elem.end(); it++){
        point_plane p;
        int x = it->first;
        int y = it->second;
        p.x = x;
        p.y = y;
        pixels_cc.push_back(p);
      }
      cc_nelem.push_back(elem.size());
    }


    //Store the pixels with the new order.
    shape *sh = get_shape(index);
    for (int i = 0; i < nb_region; i++){
      int x = pixels_cc[i].x;
      int y = pixels_cc[i].y;
      sh->pixels[i].x = x;
      sh->pixels[i].y = y;
    }

    int cc_n = cc_nelem.size();
    sh->cc_n = cc_n;
    sh->cc_of = (of_value *) malloc(cc_n*sizeof(of_value));
    sh->cc_ini = (point_plane **) malloc(cc_n*sizeof(point_plane *));
    sh->cc_nelem = (int *) malloc(cc_n*sizeof(int));


    int aux = 0;
    // printf("CC:%d\n", cc_nelem.size());
    //Crea los punteros que indican donde empieza cada componente conexa.
    for (std::vector<int>::size_type i = 0; i != cc_nelem.size(); i++){
      if (cc_nelem[i]> 1){
        total_CC_mayores++;
        // printf("CC_elemen:%d\n",cc_nelem[i] );
      }
      total_CC++;
      sh->cc_ini[i] = sh->pixels + aux;
      aux += cc_nelem[i];
      sh->cc_nelem[i] = cc_nelem[i];
    }
    // printf("Reordenados privates pixels\n");
  }//for (int index = 0; index < nb_shapes; index++)
  //TODO: Create a define with Metrics
  // printf("Porcenatje CC (pixels>1):%f\n",(1.0*total_CC_mayores/total_CC));
}

//Determine the topological  neigbors (in the cc sense) for each cc.
//TODO: Borrar printf de comprombacion
void FLSTUtility::create_topological_neigbors_list(){

  int nb_shapes = get_number_of_shapes();
  int w = _pTree->ncol;
  int h = _pTree->nrow;

  int neighborhood[8][2] = {
    {0,1}, {0,-1},{1,0},{-1,0},
    {1,1}, {1,-1},{-1,1}, {-1,-1}};
  int n_neigh = 8;

  for (int index = 0; index < nb_shapes; index++){
    const int cc_n = get_number_of_cc(index);
    shape *sh = get_shape(index);

    sh->tp_list = (tp_neighbor **)malloc(cc_n*sizeof(tp_neighbor *));
    sh->tp_n   = (int *)malloc(cc_n*sizeof(int));
    //For each connected component obtain its topologial cc.
    for (int l = 0; l < cc_n; l++){

      std::set<struct point_plane, point_comparator> pixel_neigh;
      std::set<struct point_plane, point_comparator> pixel_cc;
      std::set<struct point_plane, point_comparator>::iterator it;

      //Obtain the pixel neighbors (x,y)
      int cc_nelem = get_number_elements_cc(index,l);
      // printf("Index, cc : (%d, %d)  Numero_elem:%d \n",index, l , cc_nelem);
      for (int j = 0; j < cc_nelem; j++){
        point_plane p;
        const int x = get_tp_coord_x_cc(index, l, j);
        const int y = get_tp_coord_y_cc(index, l, j);

        for (int k =0; k < n_neigh; k++){
          int px = x + neighborhood[k][0];
          int py = y + neighborhood[k][1];
          if (px >= 0 && px < w && py >=0 && py < h){
            p.x = px;
            p.y = py;
            pixel_neigh.insert(p);
          }
        }
        p.x = x;
        p.y = y;
        pixel_cc.insert(p);
      }
      //Delete posible pixel that belongs to the cc.
      std::set<struct point_plane, point_comparator> diff;
      std::set_difference(pixel_neigh.begin(), pixel_neigh.end(), 
                          pixel_cc.begin(), pixel_cc.end(),
                          std::inserter(diff,diff.begin()), point_comparator());


      //Obtain a list that contains (index, cc) for each pixel. Only contains (index, cc) once.
      std::set<struct point_plane, point_comparator> cc_neigh;
      for (it = diff.begin(); it!=diff.end(); ++it){
        const int x = it->x;
        const int y = it->y;
        const int index = _cc_pixel_map[y*w + x].first;
        const int cc = _cc_pixel_map[y*w + x].second;
        if (index<0){
          printf("Al crear la lista interna\n");
        }
        if (cc<0){
          printf("Al crear la lista interna\n");
        }
        point_plane p;
        p.x = index;
        p.y = cc;

        if (p.x<0){
          printf("Point_x:Al crear la lista interna :%d, %d,\n",index, p.x);
        }
        if (p.y<0){
          printf("Point_y:Al crear la lista interna\n");
        }
        cc_neigh.insert(p);
      }

      // //Copy that list to each shape.
      sh->tp_n[l] = (int) cc_neigh.size();
      sh->tp_list[l] = (tp_neighbor *)malloc(cc_neigh.size()*sizeof(tp_neighbor));
      int k = 0;
      for (it = cc_neigh.begin(); it!=cc_neigh.end(); ++it){
        const int index = it->x;
        const int cc = it->y;
        if (index<0){
          printf("Index: Al recorrer cc_neigh\n");
        }
        if (cc<0){
          printf("CC: Al recorrer cc_neigh\n");
        }
        sh->tp_list[l][k].index = index;
        sh->tp_list[l][k].cc = cc;
        k++;
      }
    }
  }
}


//Get the (index,cc) related to a pixel 
void FLSTUtility::create_cc_pixel_map(){
  int nb_shapes = get_number_of_shapes();
  int w = _pTree->ncol;
  int h = _pTree->nrow;

  _cc_pixel_map.resize(w*h);
  for (int index = 0; index < nb_shapes; index++){
    const int cc_n = get_number_of_cc(index);
    //For each connected component obtain its topologial cc.
    for (int l = 0; l < cc_n; l++){
      //Obtain the pixel neighbors (x,y)
      int cc_nelem = get_number_elements_cc(index,l);
      std::pair<int,int> elem(index,l);
      if ((index < 0) || (l < 0))
        printf("_cc_pixel_map: Mierda\n");
      for (int j = 0; j < cc_nelem; j++){
        const int x = get_tp_coord_x_cc(index, l, j);
        const int y = get_tp_coord_y_cc(index, l, j);
        _cc_pixel_map[y*w + x] = elem;
      }
    }
  }
}

//Create mininum bounding box to contain each cc.
void FLSTUtility::create_min_bounding_box_cc(){
  int nb_shapes = get_number_of_shapes();
  int w = _pTree->ncol;
  int h = _pTree->nrow;

  for (int index = 0; index < nb_shapes; index++){
    const int cc_n = get_number_of_cc(index);
    shape *sh = get_shape(index);
 
    sh->x_min = (int *)malloc(cc_n*sizeof(int));
    sh->y_min = (int *)malloc(cc_n*sizeof(int));
    sh->x_max = (int *)malloc(cc_n*sizeof(int));
    sh->y_max = (int *)malloc(cc_n*sizeof(int));

    //For each cc obtain the minimum square patch that contain it.
    for (int l = 0; l < cc_n; l++){
      int cc_nelem = get_number_elements_cc(index,l);
      sh->x_min[l] = w;
      sh->y_min[l] = h;
      sh->x_max[l] = -1;
      sh->y_max[l] = -1;

      for (int j = 0; j < cc_nelem; j++){
        const int x = get_tp_coord_x_cc(index, l, j);
        const int y = get_tp_coord_y_cc(index, l, j);

        if (x <= sh->x_min[l])
            sh->x_min[l] = x;
        
        if (y <= sh->y_min[l])
            sh->y_min[l] = y;
        
        if (x >= sh->x_max[l])
            sh->x_max[l] = x;
        
        if (y >= sh->y_max[l])
            sh->y_max[l] = y;
      }
    }
  }// for(int index = 0;....)
    printf("Min BB\n");
}

//Create extended bounding box based on a radius (r).
void FLSTUtility::create_extended_bounding_box_cc(int r){
  int nb_shapes = get_number_of_shapes();
  int w = _pTree->ncol;
  int h = _pTree->nrow;

  for (int index = 0; index < nb_shapes; index++){
    const int cc_n = get_number_of_cc(index);
    shape *sh = get_shape(index);

    //For each connected component expands its boundaries.
    for (int l = 0; l < cc_n; l++){
      sh->x_min[l] = ((sh->x_min[l] - r) >=  0) ? (sh->x_min[l] - r): 0;
      sh->y_min[l] = ((sh->y_min[l] - r) >=  0) ? (sh->y_min[l] - r): 0;
      sh->x_max[l] = ((sh->x_max[l] + r) <=  w) ? (sh->x_max[l] + r + 1): w;
      sh->y_max[l] = ((sh->y_max[l] + r) <=  h) ? (sh->y_max[l] + r + 1): h;
    }
  }// for(int index = 0;....)
    printf("Max BB\n");
}

void FLSTUtility::create_extended_bounding_box_cc(int r, int index, int cc){
  int w = _pTree->ncol;
  int h = _pTree->nrow;

  shape *sh = get_shape(index);
  sh->x_min[cc] = ((sh->x_min[cc] - r) >=  0) ? (sh->x_min[cc] - r): 0;
  sh->y_min[cc] = ((sh->y_min[cc] - r) >=  0) ? (sh->y_min[cc] - r): 0;
  sh->x_max[cc] = ((sh->x_max[cc] + r) <=  w) ? (sh->x_max[cc] + r ): w;
  sh->y_max[cc] = ((sh->y_max[cc] + r) <=  h) ? (sh->y_max[cc] + r ): h;

}

//TODO: Rehacer los Makefile e incluir esto.
//TODO: REvisar como es la imagen, si ya esta en lab space.
//TODO: sigma_x, sigma_y  deber ser un valor dependiente del tamaño del patch.
//Create extended bounding box based on a radius.
void FLSTUtility::create_gaussian_weigh_bounding_box(FixedImage<float> source,
                                                     float sigma_c){
  int nb_shapes = get_number_of_shapes();
  int channels = source.get_number_of_channels();
  int w = source.get_size_x();
  int h = source.get_size_y(); 


   printf("C:%d W:%d H:%d \n",channels, w, h);
  //TODO: Create a gaussian weight for the cc == 1
  float sigma_1 = std::sqrt((float) _r_bb + 0.5);
  sigma_1 = 1.0;
  Image<float> gaussian_1 = GaussianWeights::calculate(1, 1, sigma_1, sigma_1);

  for (int index = 0; index < nb_shapes; index++){
    const int cc_n = get_number_of_cc(index);
    // printf("Total:%d, val:%d CC_N:%d\n",nb_shapes, index, cc_n);
    shape *sh = get_shape(index);
    sh->gaussian_weight = (float **)malloc(cc_n*sizeof(float *));

    //TODO: Direct access
    //For each connected component obtain, obtain the gaussian weights.
    for (int l = 0; l < cc_n; l++){
      int cc_nelem = get_number_elements_cc(index,l);
      std::vector<float> mean_val(channels,0.0);
      const int w_cc = sh->x_max[l] - sh->x_min[l];
      const int h_cc = sh->y_max[l] - sh->y_min[l];
      // printf("W: %d H: %d\n",w_cc, h_cc);
      const int size_cc = w_cc*h_cc;

      // printf("  NUm:%d, CC_elemen:%d W:%d, H:%d\n", l, cc_nelem, w_cc, h_cc);

      //Get the average value of cc for each channel 
      for (int j = 0; j < cc_nelem; j++){
        const int x = get_tp_coord_x_cc(index, l, j);
        const int y = get_tp_coord_y_cc(index, l, j);

        for (int k = 0; k < channels; k++){
          mean_val[k] += source(x, y, k);
        }
      }
      for (int k = 0; k < channels; k++){
          mean_val[k] /= (cc_nelem*1.0);
      }

      float sigma_y = (float) h_cc/2;
      float sigma_x = (float) w_cc/2;

      Image<float> gaussian;
      if ((w_cc == 1) && (h_cc == 1)){
        gaussian = gaussian_1;
      }else{
        gaussian = GaussianWeights::calculate_no_normalized(w_cc, h_cc, sigma_x, sigma_y);
      }
      sh->gaussian_weight[l] = (float *)malloc(size_cc*sizeof(float));
      //Color Gaussian
      for (int j = sh->y_min[l]; j < sh->y_max[l]; j++){
      for (int i = sh->x_min[l]; i < sh->x_max[l]; i++){
        float color = 0.0;
        for (int k = 0; k < channels; k++){
          color += pow((source(i,j,k)-mean_val[k]),2);
        }
        gaussian(i-sh->x_min[l],j-sh->y_min[l]) *= exp(-color / (2 * pow(sigma_c, 2)));
        //Copy in a non efficient way.
        sh->gaussian_weight[l][(j-sh->y_min[l])*w_cc + (i-sh->x_min[l])] = 
                                        gaussian(i-sh->x_min[l],j-sh->y_min[l]);
      }  
      }   
      //Put to one the elements which belong to the cc.
      for (int j = 0; j < cc_nelem; j++){
        const int x = get_tp_coord_x_cc(index, l, j);
        const int y = get_tp_coord_y_cc(index, l, j);
        sh->gaussian_weight[l][(y-sh->y_min[l])*w_cc + (x-sh->x_min[l])] = 1.0;
      }

    } // for (int l = ; l < cc_n; ...)
  }// for(int index = 0; index < nb_shapes;...)
}


void FLSTUtility::create_gaussian_weigh_bounding_box_bilateral(FixedImage<float> source,
                                                  float sigma_c){
    Image<float> filtered = BilateralFilter::calculate(source, 2, 1, 1, 10,  1);
    create_gaussian_weigh_bounding_box(filtered, sigma_c);
    printf("Pesos gaussianos creados\n");

}

void FLSTUtility::liberate_bounding_box(){
  int nb_shapes = get_number_of_shapes();
  for (int index = 0; index < nb_shapes; index++){
    const int cc_n = get_number_of_cc(index);
    shape *sh = get_shape(index);
    for (int l = 0; l < cc_n; l++){
      free(sh->gaussian_weight[l]);
    }
    free(sh->gaussian_weight);

    free(sh->x_min);
    free(sh->y_min);
    free(sh->x_max);
    free(sh->y_max);
  }
}

void FLSTUtility::create_bilateral_image(FixedImage<float> source){


  Image<float> filtered = BilateralFilter::calculate(source, 2, 1, 1, 10,  1);
  // filtered = IOUtility::lab_to_rgb(filtered);

  IOUtility::write_rgb_image("filtered_" + std::to_string(_pMinArea) +  ".png", filtered);
}


// ////////////////////////////Coloma CODE////////////////////////////////////////
void FLSTUtility::barOmega(int xO, int yO, int etiquetaCC,
                        double hole_graylev, int  *original,
                        int *Omegah, int dx,int dy)
{
  unsigned char eigth;

  /* LET'S COMPUTE THE HOLE Omega */
  // eigth=0; // 4-conec
  eigth = 1; // 8-conec
  hole(xO, yO, etiquetaCC, hole_graylev, hole_graylev, original, Omegah, eigth, dx,dy);
}


unsigned char FLSTUtility::hole(int xp, int yp, int etiquetaCC, 
                  double ming, double maxg,int  *Image, int *OmegaH, 
                  unsigned char connecty8, int dx, int dy)
{
  unsigned char following;
  long p;

  p=(long)dx; p*=(long)yp; p+=(long)xp;
  OmegaH[p]=etiquetaCC;

  following=1;

  if (((xp-1) >= 0))
  { 
    if ((OmegaH[p-1] == 0) && (Image[p-1] >= ming) && (Image[p-1] <= maxg))
    {
      following=hole(xp-1, yp, etiquetaCC, ming, maxg, Image, OmegaH, connecty8, dx,dy);
    }
  }

  if ((following) && ((yp-1) >= 0))
  {
    if ((OmegaH[p-dx] == 0) && (Image[p-dx] >= ming) && (Image[p-dx] <= maxg))
    {
      following=hole(xp, yp-1, etiquetaCC, ming, maxg, Image, OmegaH, connecty8, dx,dy);
    }
  }

  if ((following) && ((xp+1) < dx))
  {
    if ((OmegaH[p+1] == 0) && (Image[p+1] >= ming) && (Image[p+1] <= maxg))
    {
      following=hole(xp+1, yp, etiquetaCC, ming, maxg, Image, OmegaH, connecty8, dx,dy);
    }
  }

  if ((following) && ((yp+1) < dy))
  {
    if ((OmegaH[p+dx] == 0) && (Image[p+dx] >= ming) && (Image[p+dx] <= maxg))
    {
      
      following=hole(xp, yp+1, etiquetaCC, ming, maxg, Image, OmegaH, connecty8, dx,dy);
      
    }
  }

 if(connecty8 == 1)
 {
  if ((following) && ((xp-1) >= 0) && ((yp-1) >= 0))
  {
    if ((OmegaH[p-1-dx] == 0) && (Image[p-1-dx] >= ming) && (Image[p-1-dx] <= maxg))
     {
       following=hole(xp-1, yp-1, etiquetaCC, ming, maxg, Image, OmegaH, connecty8, dx,dy);
     }
  }

  if ((following) && ((yp-1) >= 0) && ((xp+1) < dx))
  {
    if ((OmegaH[p-dx+1] == 0) && (Image[p-dx+1] >= ming) && (Image[p-dx+1] <= maxg))
    {
      following=hole(xp+1, yp-1, etiquetaCC, ming, maxg, Image, OmegaH, connecty8, dx, dy);
    }
  }

  if ((following) && ((xp+1) < dx) && ((yp+1) < dy))
  {
    if ((OmegaH[p+1+dx] == 0) && (Image[p+1+dx] >= ming) && (Image[p+1+dx] <= maxg))
    {
      following=hole(xp+1, yp+1, etiquetaCC, ming, maxg, Image, OmegaH, connecty8, dx,dy);
    }
  }

  if ((following) && ((yp+1) < dy) && ((xp-1) >= 0))
  {
    if ((OmegaH[p+dx-1] == 0) && (Image[p+dx-1] >= ming) && (Image[p+dx-1] <= maxg))
    {
      following=hole(xp-1, yp+1, etiquetaCC, ming, maxg, Image, OmegaH, connecty8, dx,dy);
    }
  }
 }

  if (following) return(1);
  else return(0);
}

////////////////////////////////////////////////////////////////////////////////

Fimage FLSTUtility::FixedImage_to_Fimage(FixedImage<float> source){
    // Image<float> input = IOUtility::lab_to_rgb(source);
    //TODO: It assumes that image format is lab space.
    Image<float> input = source;
    const float *im = input.raw();
    int channels = source.get_number_of_channels();  
    printf("number of channel of source : %d \n", channels); /*always 3*/

    uint w = source.get_size_x();
    uint h = source.get_size_y();
    Fimage my_image = mw_new_fimage();
    mw_change_fimage(my_image, h, w);

    std::vector<float> image_tmp(w*h);

    for (uint y = 0; y < h; y++) {
    for (uint x = 0; x < w; x++) {
      int index = channels * (y * w + x);
      //TODO: It assumes that image format is lab space. We only take the
      // luminance channel.
      // if (channels == 3){
      //    image_tmp[y*w +x] = .2126*im[index] + .7152*im[index + 1] + .0722*im[index +2];
      // }
      if (channels == 3){
         image_tmp[y*w +x] = im[index];  
      }
      else{
        image_tmp[y*w +x]  = im[index];   
      }
    }
    }

    image_normalization(image_tmp.data(),image_tmp.data(),w*h);
    for (uint y = 0; y < h; y++) {
    for (uint x = 0; x < w; x++) {
       my_image->gray[y*w + x] = image_tmp[y*w +x];
    }
    }


    return my_image;  
}

FixedImage<float> FLSTUtility::Fimage_to_FixedImage(Fimage source){

  uint w = source->ncol;
  uint h = source->nrow;
  float *im = source->gray;


  Image<float> image(w, h);  /*1 channel*/

  // copy from image_data to Image<float>
  for (uint y = 0; y < h; y++) {
    for (uint x = 0; x < w; x++) {
      image(x,y) = im[y*w + x];
    }
  }

  return image;  
}


/*Test for DEBUG*/
Image<float> FLSTUtility::test_library(){

    Fimage im_out = mw_new_fimage();

    flst_reconstruct(_pTree, im_out);
    prune_cc(im_out) ; 
    Image<float> output = Fimage_to_FixedImage(im_out);
    mw_delete_fimage(im_out);
    return output;
}

Image<float> FLSTUtility::test_library(FixedImage<float> source, int pMinArea){

    Fimage input = FixedImage_to_Fimage(source);
    Fimage im_out = mw_new_fimage();
    Shapes pTree = mw_new_shapes();

    flst(&pMinArea, input, pTree);
    flst_reconstruct(pTree, im_out);
    prune_cc(im_out);
    Image<float> output = Fimage_to_FixedImage(im_out);

    mw_delete_shapes(pTree);
    mw_delete_fimage(im_out);
    mw_delete_fimage(input);

    return output;
}



//Draws the area, the region and the cc of a certain shape in three different png
void FLSTUtility::draw_cc_shape(int index){

  int w = _pTree->ncol;
  int h = _pTree->nrow;


  Image<float> regiones(w,h,(uint)3);

  int color_n = 12;
  std::vector<float> color_list(color_n*3);
    //Red
  color_list[0] = 255;
  color_list[1] = 0;
  color_list[2] = 0;
  //Green
  color_list[3] = 0;
  color_list[4] = 255;
  color_list[5] = 0;
  //Blue
  color_list[6] = 0;
  color_list[7] = 0;
  color_list[8] = 255;
  //Other
  color_list[9]  = 0;
  color_list[10] = 255;
  color_list[11] = 255;
  //Other
  color_list[12] = 255;
  color_list[13] = 255;
  color_list[14] = 0;
  //Other
  color_list[15] = 255;
  color_list[16] = 128;
  color_list[17] = 0;
  //Other
  color_list[18] = 102;
  color_list[19] = 102;
  color_list[20] = 0;
  //Other
  color_list[21] = 127;
  color_list[22] = 0;
  color_list[23] = 255;
  //Other
  color_list[24] = 255;
  color_list[25] = 0;
  color_list[26] = 127;
  //Morado
  color_list[27] = 255;
  color_list[28] = 0;
  color_list[29] = 255;
  //WHite
  color_list[30] = 255;
  color_list[31] = 255;
  color_list[32] = 255;
  //Gray
  color_list[33] = 128;
  color_list[34] = 128;
  color_list[35] = 128;

  //Background in black
  for (int y = 0; y < h; y++)
  for (int x = 0; x < w; x++){
    regiones(x,y, 0) = 0.0;
    regiones(x,y, 1) = 0.0;
    regiones(x,y, 2) = 0.0;
  }
  printf("Index: %d\n", index);
  //Area gray
  int area = get_area(index);
  printf("  Area: %d\n", area);
  for (int l = 0; l < area; l++){
    const int x = _root_to_leave[index]->pixels[l].x;
    const int y = _root_to_leave[index]->pixels[l].y;
    // printf("     Element: %d  (%d,%d)\n", l,x,y);
    regiones(x,y, 0) = color_list[33];
    regiones(x,y, 1) = color_list[34];
    regiones(x,y, 2) = color_list[35];
  }
  IOUtility::write_rgb_image("Shape_" + std::to_string(index) + ".png", regiones);
  
  //Region in White
  int region = get_region(index);
  printf("  Region: %d\n", region);
  for (int l = 0; l < region; l++){
    const int x = get_tp_coord_x(index, l);
    const int y = get_tp_coord_y(index, l);
    regiones(x,y, 0) = color_list[30];
    regiones(x,y, 1) = color_list[31];
    regiones(x,y, 2) = color_list[32];
  }
  IOUtility::write_rgb_image("Region_" + std::to_string(index) + ".png", regiones);
  

  //CC in color
  int cc_n = get_number_of_cc(index);
  printf("CC_N:%d\n", cc_n);
  for (int l = 0; l < cc_n; l++){
    int cc_nelem = get_number_elements_cc(index, l);
      printf("   inner_elem:%d\n", cc_nelem);
    for (int j = 0; j < cc_nelem; j++){
        const int x = get_tp_coord_x_cc(index, l, j);
        const int y = get_tp_coord_y_cc(index, l, j);
        if (l < 9){
          regiones(x,y,0) = color_list[l*3 + 0];
          regiones(x,y,1) = color_list[l*3 + 1];
          regiones(x,y,2) = color_list[l*3 + 2];
        }else{
          regiones(x,y,0) = color_list[27];
          regiones(x,y,1) = color_list[28];
          regiones(x,y,2) = color_list[29];
        }
    }
  }
  IOUtility::write_rgb_image("CC_" + std::to_string(index) + ".png", regiones);
}



//Draws the area, the region and the cc of a certain shape in three different png
void FLSTUtility::draw_topological_neighbors_cc(int index){

  int w = _pTree->ncol;
  int h = _pTree->nrow;


  Image<float> regiones(w,h,(uint)3);

  int color_n = 3;
  std::vector<float> color_list(color_n*3);
  //Red
  color_list[0] = 255;
  color_list[1] = 0;
  color_list[2] = 0;
  //Morado
  color_list[3] = 255;
  color_list[4] = 0;
  color_list[5] = 255;
  //WHite
  color_list[6] = 255;
  color_list[7] = 255;
  color_list[8] = 255;
  //Gray
  color_list[9] = 128;
  color_list[10] = 128;
  color_list[11] = 128;

  //Background in black
  for (int y = 0; y < h; y++)
  for (int x = 0; x < w; x++){
    regiones(x,y, 0) = 0.0;
    regiones(x,y, 1) = 0.0;
    regiones(x,y, 2) = 0.0;
  }
  printf("Index: %d\n", index);
  //Area gray
  int area = get_area(index);
  printf("  Area: %d\n", area);
  for (int l = 0; l < area; l++){
    const int x = _root_to_leave[index]->pixels[l].x;
    const int y = _root_to_leave[index]->pixels[l].y;
    // printf("     Element: %d  (%d,%d)\n", l,x,y);
    regiones(x,y, 0) = color_list[9];
    regiones(x,y, 1) = color_list[10];
    regiones(x,y, 2) = color_list[11];
  }

  //Region in White
  int region = get_region(index);
  printf("  Region: %d\n", region);
  for (int l = 0; l < region; l++){
    const int x = get_tp_coord_x(index, l);
    const int y = get_tp_coord_y(index, l);
    regiones(x,y, 0) = color_list[6];
    regiones(x,y, 1) = color_list[7];
    regiones(x,y, 2) = color_list[8];
  }
  IOUtility::write_rgb_image("Region_" + std::to_string(index) + ".png", regiones);
  

  int cc_n = get_number_of_cc(index);
  printf("CC_N:%d\n", cc_n);
  for (int l = 0; l < cc_n; l++){
    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++){
      regiones(x,y, 0) = 0.0;
      regiones(x,y, 1) = 0.0;
      regiones(x,y, 2) = 0.0;
    }
    //Region in White
    int region = get_region(index);
    printf("  Region: %d\n", region);
    for (int l = 0; l < region; l++){
      const int x = get_tp_coord_x(index, l);
      const int y = get_tp_coord_y(index, l);
      regiones(x,y, 0) = color_list[6];
      regiones(x,y, 1) = color_list[7];
      regiones(x,y, 2) = color_list[8];
    }

    //Neighbors in color
    int cc_nelem = get_number_elements_cc(index, l);
      printf("   inner_elem:%d\n", cc_nelem);
    for (int j = 0; j < cc_nelem; j++){
      const int x = get_tp_coord_x_cc(index, l, j);
      const int y = get_tp_coord_y_cc(index, l, j);
        regiones(x,y,0) = color_list[0];
        regiones(x,y,1) = color_list[1];
        regiones(x,y,2) = color_list[2];
    }

    shape *sh = get_shape(index);
    int neig_n = sh->tp_n[l];
    for (int i = 0; i < neig_n; i++){
      const int neig_index = sh->tp_list[l][i].index;
      const int neig_cc = sh->tp_list[l][i].cc;
      int cc_nelem = get_number_elements_cc(neig_index, neig_cc);
      printf("Topological Neig::%d\n", cc_nelem);
      int color_rand[3]  = {rand() %255, rand() %255, rand() %255};
      for (int j = 0; j < cc_nelem; j++){
        const int x = get_tp_coord_x_cc(neig_index, neig_cc, j);
        const int y = get_tp_coord_y_cc(neig_index, neig_cc, j);
        regiones(x,y,0) = color_rand[0];
        regiones(x,y,1) = color_rand[1];
        regiones(x,y,2) = color_rand[2];
      }
    }
    IOUtility::write_rgb_image("NeighborCC_" + std::to_string(index) +  "_" + std::to_string(l) + ".png", regiones);
  }
}


//Draws a png where each conected component is filled with a different color.
void FLSTUtility::draw_all_cc(){
  int nb_shapes = get_number_of_shapes();
  int w = _pTree->ncol;
  int h = _pTree->nrow;

  Image<float> regiones(w,h,(uint)3);  /*3 channels*/
  printf("number of channels of regiones : %d \n" , regiones.get_number_of_channels()) ; /*Lucie*/


  for (int index = 0; index < nb_shapes; index++){
    const int cc_n = get_number_of_cc(index);
    //For each cc obtain the minimum square patch that contain it.
    for (int l = 0; l < cc_n; l++){
      int cc_nelem = get_number_elements_cc(index,l);
      int color_rand[3]  = {rand() %255, rand() %255, rand() %255};

      for (int j = 0; j < cc_nelem; j++){
        const int x = get_tp_coord_x_cc(index, l, j);
        const int y = get_tp_coord_y_cc(index, l, j);
        regiones(x,y,0) = color_rand[0];
        regiones(x,y,1) = color_rand[1];
        regiones(x,y,2) = color_rand[2];
      }
    }
  }// for(int index = 0;....)



   IOUtility::write_rgb_image("AllCC_" + std::to_string(_pMinArea) +  ".png", regiones);
}


//Check how many pixels have an area below minArea
void FLSTUtility::test_region_pixels(int minArea){
  int nb_shapes = _pTree->nb_shapes;
  int w = _pTree->ncol;
  int h = _pTree->nrow;
  int j = 0;
  for (int i = 0; i < nb_shapes; i++){
      const int region = _root_to_leave[i]->region;
      if (region < minArea)
          j++;
  }
  printf("Porcentaje:%f\n",(1.0*j)/(w*h) );
}


//Check that each shape has at least one private pixel.
void FLSTUtility::test_at_least_one_pixel()
{
  shape *sh;

  for (int k = 0; k < _pTree->nb_shapes; k++){
    sh = &_pTree->the_shapes[k];
    shape *child;
    int n_pixels_childs = 0;
    for (child = sh->child; child!=NULL; child = child->next_sibling){     
      n_pixels_childs += child->area;
    }
    if (sh->area == n_pixels_childs){
      printf("Menor Area y REgion iguales\n");
    }
  }

}

//Check that te process to obtain the pixel list in a efficient way is equivalent to do the non-efficient way.
void FLSTUtility::test_shape_without_holes(){

  int k = 0;
  int diferentes = 0;
  int iguales = 0;
  shape *sh;
  for (k = 0; k < _pTree->nb_shapes; k++){
    // printf("Total:%d, val:%d       ",_pTree->nb_shapes,k );
    //Version eficiente
    sh = &(_pTree->the_shapes[k]);
    point_plane *child_pixels;
    if (sh->child!=NULL){
      child_pixels = sh->child->pixels;
      point_plane *aux;
      int j = 0;
      aux = sh->pixels;
      while (aux!=child_pixels){
        aux++;
        j++;
      }
      sh->region = j;

     // //Version No-eficiente
        int n_pixels_childs = 0;
        int n_pixels_shape = sh->area;
            
        std::vector<struct point_plane> v_origin(n_pixels_shape);
        std::vector<struct point_plane> diff(n_pixels_shape);
        std::vector<struct point_plane>::iterator it;
       
        int i = 0;
        for (i= 0;i < n_pixels_shape; i++){
          v_origin[i].x = sh->pixels[i].x;
          v_origin[i].y = sh->pixels[i].y;
        }

        shape *child;
        n_pixels_childs = 0;
        for (child = sh->child; child!=NULL; child = child->next_sibling){
        // for (child = sh->child; child!=NULL; child = child->next_sibling){     
          n_pixels_childs += child->area;
        }
        std::vector<struct point_plane> v_child(n_pixels_childs);   
 
        int m = 0;
        for (child = sh->child; child!=NULL; child = child->next_sibling){     
          for (int l = 0;l < child->area; l++, m++){
            v_child[m].x = child->pixels[l].x;
            v_child[m].y = child->pixels[l].y;
          }
        }
        if (m!=n_pixels_childs)
          printf(" \nCosa m:%d n_pixels_childs:%d\n", m, n_pixels_childs);

        std::sort(v_origin.begin(), v_origin.begin() + n_pixels_shape, point_comparator());
        std::sort(v_child.begin(),  v_child.begin() + n_pixels_childs, point_comparator());


        it = std::set_difference(v_origin.begin(), v_origin.end(), 
                                 v_child.begin(),   v_child.end(),
                                diff.begin(), point_comparator());

        diff.resize(it-diff.begin());
               
        //Comprobacion son lo mismo
        std::vector<struct point_plane> effi(sh->region);
        for (i = 0; i < sh->region; i++){
          effi[i].x = sh->pixels[i].x;
          effi[i].y = sh->pixels[i].y;
        }
        std::sort(effi.begin(),effi.end(), point_comparator());
        std::sort(diff.begin(),diff.end(), point_comparator());
        bool flag =
          std::equal(effi.begin(), effi.end(), 
                     diff.begin(), equal_comparator());
        if (!flag){
          diferentes++;
          printf("No son iguales $$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
          printf("Origin\n");        
          for (it = v_origin.begin(); it!=v_origin.end(); it++){
            printf(" (%d,%d) ", it->x,it->y);
          }
          printf("\n");
          printf("Child\n");        
          for (it = v_child.begin(); it!=v_child.end(); it++){
            printf(" (%d,%d) ", it->x,it->y);
          }
          printf("\n");
        }else{
            iguales++;
          }
    }else{
      iguales++;
      sh->region = sh->area;
    }
  }
  printf("Total Holes: %f Distintos HOles %f \n",(1.0*iguales)/_pTree->nb_shapes, (1.0*diferentes)/_pTree->nb_shapes);
}

//Check if all there is no missing pixel when you obtain the pixel list using regions.
void FLSTUtility::test_no_repetition()
{
  shape *sh;
  int w = _pTree->ncol;
  int h = _pTree->nrow;
  int k;
  printf("Total:%d\n",w*h);
  std::vector<struct point_plane> v_origin;
  std::vector<struct point_plane> v_shape;
  for (int j = 0; j < h; j++){
    for (int i = 0; i < w; i++){
      point_plane p;
      p.x = i;
      p.y = j;
      v_origin.push_back(p);
    }
  }

  for (k = 0; k < _pTree->nb_shapes; k++){
    sh = &_pTree->the_shapes[k];
    for (int l = 0; l < sh->region; l++){
      v_shape.push_back(sh->pixels[l]);
    }
  }

  if (v_origin.size()!= v_shape.size()){
    printf("No tienen el mismo tamaño v_origin:%lu y v_shape:%lu\n",v_origin.size(), v_shape.size());
  }

  std::sort(v_origin.begin(), v_origin.end(), point_comparator());
  std::sort(v_shape.begin(),  v_shape.end(),  point_comparator());
  bool flag = std::equal(v_origin.begin(), v_origin.end(), 
                     v_shape.begin(), equal_comparator());

  if (flag){
    printf("Son los mismos elementos\n");
  }else{
    printf("No hago bien lo de las regiones\n");
  }
}



void FLSTUtility::getminmax(
  float *min,     // output min
  float *max,     // output max
  float *x, // input array
  int n           // array size
)
{
  *min = *max = x[0];
  for (int i = 1; i < n; i++) {
    if (x[i] < *min)
      *min = x[i];
    if (x[i] > *max)
      *max = x[i];
  }
}
/**
 *
 * Function to normalize the images between 0 and 255
 *
 **/
void FLSTUtility::image_normalization(
    float *a,  // input image0
    float *an,       // normalized output image0
    int size          // size of the image
    )
{
  float max, min;
  // obtain the max and min of both images
  getminmax(&min, &max, a, size);
  const float den = max - min;

  if (den > 0)
    // normalize both images
    for (int i = 0; i < size; i++)
    {
      an[i] = std::floor((255*(a[i] - min) / den));
    }

  else
    // copy the original images
    for (int i = 0; i < size; i++)
    {
      an[i] = a[i];
    }
}



//Total Variation, global and local, for an auxiliary image  /*Lucie*/

float FLSTUtility::global_variation(Fimage auxImage)
 {
  int ncol = auxImage->ncol ;
  int nrow = auxImage->nrow ;
  int x ; 
  float gradx, grady, TV = 0 ; 

  for (int i = 1 ; i < nrow - 1 ; i++){
    for (int j = 1 ; j < ncol - 1 ; j++){
      x = i*ncol + j ;
      gradx = (auxImage->gray[x+1] - auxImage->gray[x-1])/2; 
      grady = (auxImage->gray[x+ncol] - auxImage->gray[x-ncol])/2 ;
      TV += abs(gradx) + abs(grady) ;
    }
  }
  return TV ; 
} 


Fimage FLSTUtility::create_patch_parent(Fimage im_out, int index, int cc)  // problem with TV 
{
  //printf("create aux Image parent start \n");
  call_parent++ ;
  int w = _pTree->ncol;
  int numPixel , numPixel_im_out ;

  create_extended_bounding_box_cc(2 , index, cc);

  int x_min = get_x_min_cc(index, cc) ;
  //printf ("x_min %d call_parent:%d \n" , x_min , call_parent );
  int x_max = get_x_max_cc(index, cc) ; 
  //printf ("x_max %d call_parent :%d \n" , x_max , call_parent );
  int y_min = get_y_min_cc(index, cc) ; 
  //printf ("y_min %d call_parent :%d \n" , y_min , call_parent );
  int y_max = get_y_max_cc(index, cc) ;
  //printf ("y_max %d call_parent :%d \n" , y_max , call_parent );

  int ncol = x_max - x_min + 1 ; 
  int nrow = y_max - y_min + 1 ;

  Fimage patch = mw_new_fimage() ; 
  mw_change_fimage(patch, nrow, ncol) ;


  for (int i = 0; i < ncol ; i++) {
    for (int j = 0 ; j < nrow  ; j++){
      numPixel = j*ncol + i ; //here ? 
      numPixel_im_out = (j+y_min)*w + (i + x_min) ;
      patch->gray[numPixel] = im_out->gray[numPixel_im_out]  ; 
    }
  }

  //to optimize

  int cc_nelem = get_number_elements_cc(index,cc);
  for (int j = 0; j < cc_nelem; j++){
    const int x = get_tp_coord_x_cc(index, cc, j) - x_min;
    const int y = get_tp_coord_y_cc(index, cc, j) - y_min;
    numPixel = y*ncol + x ; 
    if (_pTree->the_shapes[index].parent == NULL) { // to optimize !!
      patch->gray[numPixel] = _pTree->the_shapes[index].value ;
    }
    else {
      patch->gray[numPixel] = _pTree->the_shapes[index].parent->value ; 
    }
  }    

  //FixedImage<float> Fpatch = Fimage_to_FixedImage(patch); //delete ?
  //IOUtility::write_mono_image("patch_parent.png", Fpatch);  


  create_extended_bounding_box_cc(-2 , index, cc);
  //printf("create aux Image parent end \n");
  return patch ; 
}

Fimage FLSTUtility::create_patch_child(Fimage im_out, int index, int cc)  //problem with TV 
{
  //printf("create aux Image child start \n");
  call_child++ ;
  int w = _pTree->ncol;
  int numPixel , numPixel_im_out ; 

  create_extended_bounding_box_cc(2 , index, cc);

  int x_min = get_x_min_cc(index, cc) ;
  //printf ("x_min %d call_child :%d \n" , x_min , call_child );
  int x_max = get_x_max_cc(index, cc) ; 
  //printf ("x_max %d call_child :%d \n" , x_max , call_child );
  int y_min = get_y_min_cc(index, cc) ; 
  //printf ("y_min %d call_child :%d \n" , y_min , call_child );
  int y_max = get_y_max_cc(index, cc) ;
  //printf ("y_max %d call_child :%d \n" , y_max , call_child );

  int ncol = x_max - x_min + 1  ; 
  int nrow = y_max - y_min + 1;

  Fimage patch = mw_new_fimage() ; 
  mw_change_fimage(patch, nrow, ncol) ;


  for (int i = 0; i < ncol  ; i++) {
    for (int j = 0 ; j < nrow  ; j++){
      numPixel = j*ncol + i ; //here ? 
      numPixel_im_out = (j+y_min)*w + (i + x_min) ;
      patch->gray[numPixel] = im_out->gray[numPixel_im_out]  ; 
    }
  }



  //to optimize

  int cc_nelem = get_number_elements_cc(index,cc);
  for (int j = 0; j < cc_nelem; j++){
    const int x = get_tp_coord_x_cc(index, cc, j) - x_min;
    const int y = get_tp_coord_y_cc(index, cc, j) - y_min ;
    numPixel = y*ncol + x ; 
    if (_pTree->the_shapes[index].child == NULL) { // to optimize !!
      patch->gray[numPixel] = _pTree->the_shapes[index].value ;
    }
    else {
      patch->gray[numPixel] = _pTree->the_shapes[index].child->value ; 
    }
  }    

  //FixedImage<float> Fpatch = Fimage_to_FixedImage(patch); 
  //IOUtility::write_mono_image("patch_child.png", Fpatch);

  create_extended_bounding_box_cc(-2 , index, cc);
  //printf("create aux Image child end \n");
  return patch ; 
}


Fimage FLSTUtility::create_auxImage_parent(int index, int cc){  //create an image of size of the image and replace the value of the cc of index by parent value
  call_parent++ ;
  int nb_shapes = get_number_of_shapes() ; 
  int w = _pTree->ncol;
  int h = _pTree->nrow;
  int numPixel ; 


  Fimage auxImage = mw_new_fimage() ; 
  mw_change_fimage(auxImage, h, w) ;

  for (int shape = 0; shape < nb_shapes; shape++){
    const int cc_n = get_number_of_cc(shape); // why const ? 
    for (int l = 0; l < cc_n; l++){
      int cc_nelem = get_number_elements_cc(shape,l);
        for (int j = 0; j < cc_nelem; j++){
          const int x = get_tp_coord_x_cc(shape, l, j);
          const int y = get_tp_coord_y_cc(shape, l, j);
          numPixel = y*w + x ; 
          if (shape != index || l != cc){
            auxImage->gray[numPixel] = _pTree->the_shapes[shape].value ; 
          }
          else { 
            if (_pTree->the_shapes[shape].parent == NULL) { // to optimize !!
              auxImage->gray[numPixel] = _pTree->the_shapes[shape].value ;
            }
            else {
              auxImage->gray[numPixel] = _pTree->the_shapes[shape].parent->value ; 
            }
          }    
        }
      
    }
  }

  //FixedImage<float> FauxImage = Fimage_to_FixedImage(auxImage); 
  //IOUtility::write_mono_image("patch.png", FauxImage);      


  return auxImage ;
}


Fimage FLSTUtility::create_auxImage_child(int index, int cc) {  //create an image of size of the image and replace the value of the cc of index by child value 
  call_child++ ; 
  int nb_shapes = get_number_of_shapes() ; 
  int w = _pTree->ncol;
  int h = _pTree->nrow;
  int numPixel ; 


  Fimage auxImage = mw_new_fimage() ; 
  mw_change_fimage(auxImage, h, w) ;

  for (int shape = 0; shape < nb_shapes; shape++){
    const int cc_n = get_number_of_cc(shape); // why const ? 
    for (int l = 0; l < cc_n; l++){
      int cc_nelem = get_number_elements_cc(shape,l);
      for (int j = 0; j < cc_nelem; j++){
        const int x = get_tp_coord_x_cc(shape, l, j);
        const int y = get_tp_coord_y_cc(shape, l, j);
        numPixel = y*w + x ; 
        if (shape != index || l != cc){
          auxImage->gray[numPixel] = _pTree->the_shapes[shape].value ; 
        }
        else { 
          if (_pTree->the_shapes[shape].child == NULL) { // to optimize !!
            auxImage->gray[numPixel] = _pTree->the_shapes[shape].value ;
          }
          else {
            auxImage->gray[numPixel] = _pTree->the_shapes[shape].child->value ; 
          }
        }      
      }  
    }
  }

  //FixedImage<float> FauxImage = Fimage_to_FixedImage(auxImage); 
  //IOUtility::write_mono_image("patch.png", FauxImage);      
  return auxImage ;
}


void FLSTUtility::prune_cc(Fimage im_out){  //modify the output Fimage pointer ?? 
  printf("Call of prune_cc \n") ; 
  int nb_shapes = get_number_of_shapes() ;
  int w = _pTree->ncol ,ncol , x_min , y_min;
  int numPixel, numPixel_auxImage;
  float global_variation_child , global_variation_parent ; 


  for (int index = 0; index < nb_shapes ; index ++){
    int cc_n = get_number_of_cc(index); 
    if (cc_n > 1){
     for (int i = 0; i < cc_n; i++ ){ 
      int cc_nelem = get_number_elements_cc(index,i);
        if (cc_nelem < _pMinArea) {
          Fimage auxImage_child = create_patch_child(im_out, index, i) ;  //create at the begining ? 
          Fimage auxImage_parent = create_patch_parent(im_out, index, i) ;
          global_variation_parent = global_variation(auxImage_parent) ;
          //printf("TV parent : %f \n ", global_variation_parent);
          global_variation_child = global_variation(auxImage_child) ;
          //printf("TV child : %f \n ", global_variation_child);
          ncol = auxImage_child->ncol ; //
          x_min =  get_x_min_cc(index, i)  ; //
          //printf("xmin prune : %d \n",x_min );
          y_min =  get_y_min_cc(index, i) ; //
          for (int j = 0; j < cc_nelem; j++){
            const int x = get_tp_coord_x_cc(index, i, j) ;
            const int y = get_tp_coord_y_cc(index, i, j) ;
            numPixel = y*w + x ; 
            numPixel_auxImage = (y + 2 - y_min)*ncol + (x + 2 - x_min) ; // 
            if (global_variation_child < global_variation_parent) {  
               im_out->gray[numPixel] = auxImage_child->gray[numPixel_auxImage] ; //
            }            
            else {
              im_out->gray[numPixel] = auxImage_parent->gray[numPixel_auxImage] ; 
            }
          } 
          mw_delete_fimage(auxImage_child);
          mw_delete_fimage(auxImage_parent);
        }
      }
    }
  }
  printf("number of call of create_auxImage_parent : %d\n", call_parent );
  printf("number of call of create_auxImage_child : %d\n", call_child );
}


Fimage FLSTUtility::prune_cc() { // create an auxiliary image but diverge with a high number of shapes 

  printf("Call of prune_cc \n") ; 
  int nb_shapes = get_number_of_shapes() ;
  int w = _pTree->ncol , h = _pTree->nrow, ncol , x_min , y_min;
  int numPixel, numPixel_auxImage;
  float global_variation_child , global_variation_parent ; 

  Fimage auxImage = mw_new_fimage() ; 
  mw_change_fimage(auxImage, h, w) ;


  for (int index = 0; index < nb_shapes ; index ++){
    int cc_n = get_number_of_cc(index); 
    if (cc_n > 1){
     for (int i = 0; i < cc_n; i++ ){ 
      int cc_nelem = get_number_elements_cc(index,i);
        if (cc_nelem < _pMinArea) {
          Fimage auxImage_child = create_patch_child(auxImage, index, i) ;  //create at the begining ? 
          Fimage auxImage_parent = create_patch_parent(auxImage, index, i) ;
          global_variation_parent = global_variation(auxImage_parent) ;
          //printf("TV parent : %f \n ", global_variation_parent);
          global_variation_child = global_variation(auxImage_child) ;
          //printf("TV child : %f \n ", global_variation_child);
          ncol = auxImage_child->ncol ; //
          x_min =  get_x_min_cc(index, i)  ; //
          //printf("xmin prune : %d \n",x_min );
          y_min =  get_y_min_cc(index, i) ; //
          for (int j = 0; j < cc_nelem; j++){
            const int x = get_tp_coord_x_cc(index, i, j) ;
            const int y = get_tp_coord_y_cc(index, i, j) ;
            numPixel = y*w + x ; 
            numPixel_auxImage = (y + 2 - y_min)*ncol + (x + 2 - x_min) ; // 
            if (global_variation_child < global_variation_parent) {  
               auxImage->gray[numPixel] = auxImage_child->gray[numPixel_auxImage] ; //
            }            
            else {
              auxImage->gray[numPixel] = auxImage_parent->gray[numPixel_auxImage] ; 
            }
          } 
          mw_delete_fimage(auxImage_child);
          mw_delete_fimage(auxImage_parent);
        }
      }
    }
  }

  printf("number of call of create_auxImage_parent : %d\n", call_parent );
  printf("number of call of create_auxImage_child : %d\n", call_child );

  return auxImage ; 
}


// void FLSTUtility::prune_cc() { // modify the tree but too complicated and not reliable 

//   int nb_shapes = get_number_of_shapes() ; 
//   int w = _pTree->ncol;
//   int h = _pTree->nrow;


//   for (int index = 0; index < nb_shapes ; index ++){
//     int cc_n = get_number_of_cc(index); 
//     if (cc_n > 1){
//      for (int i = 0; i < cc_n; i++ ){ 
//       int cc_nelem = get_number_elements_cc(index,i);
//         if (cc_nelem < _pMinArea) {
//           Fimage auxImage_parent = create_auxImage_parent(index, i) ;
//           Fimage auxImage_child = create_auxImage_child(index, i) ;  
//           if (global_variation(auxImage_child) < global_variation(auxImage_parent)) {  
//             for (int j = 0; j < cc_nelem; j++){
//               if (_pTree->the_shapes[index].child == NULL) {
//                 return ; 
//               }
//               else {
//                 modify_tree() ;
//               }
//             }
//           }
//           else {
//             for (int j = 0; j < cc_nelem; j++){
//               if (_pTree->the_shapes[index].parent == NULL) {
//                 return;
//               }
//               else {
//                 modify_tree() ; 
//               }
//             }
//           } 
//         }
//       }
//     }  
//   }

//   return ; 
// }



