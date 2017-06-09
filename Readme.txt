 Si te da problemas para compilar en el make, compila con esto

 CC=clang-omp CXX=clang-omp++ cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc-5 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-5 ..


 CC=clang-omp CXX=clang-omp++ cmake ..



 #Información sobre lo que hace cada ejecutable

 main_flst - Ejecutable que presenta el mismo comportamiento que el código de flst, encapsulado y ejecutado a través de una clase.
 main_patchmatch - Ejecutable que ejecuta nnf mediante patchmatch utilizando el código de Vadim.

 main_patchmatch_tree - Nuestro código qeu intenta incorporar el código del arbol de formas.