cmake_minimum_required (VERSION 2.6) #specify the required version of cmake
project (IIO)                        #give a name to the project/solution

# ADD LOCAL PATH TO SEARCH MODULES
list (APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_PREFIX}")
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})



# SETUP A MAKE CLEAN TO CLEAN EVERYTHING
#
#
set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES 
   "iioConfig.cmake;CMakeCache.txt;CMakeFiles;Makefile;cmake_install.cmake;iio;iion;iio.o")


# LIBRARIES
#
#

#set (LIBS "")



# OTHER LIBRARIES
FIND_PACKAGE(TIFF)
IF(TIFF_FOUND)
   include_directories(${TIFF_INCLUDE_DIR})
   #  link_directories(${TIFF_LIBRARY})
   SET( LIBS ${LIBS} ${TIFF_LIBRARIES})
ENDIF()

FIND_PACKAGE(JPEG)
if(JPEG_FOUND)
   include_directories(${JPEG_INCLUDE_DIR})
   # link_directories(${JPEG_LIBRARY})
   SET( LIBS ${LIBS} ${JPEG_LIBRARIES})
ENDIF()

FIND_PACKAGE(PNG)
IF(PNG_FOUND)
   include_directories(${PNG_INCLUDE_DIR})
   #  link_directories(${PNG_LIBRARY})
   SET( LIBS ${LIBS} ${PNG_LIBRARIES})
ENDIF()


# DETECT LINUX
#
#
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   # Linux specific code flags
   SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DI_CAN_HAS_LINUX -D_GNU_SOURCE")
ENDIF()


# IIO
#
#
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -std=c99 -O3 -D_GNU_SOURCE")
message( "FLAGS: ${CMAKE_C_FLAGS}")
set(FLST_LIB_SRCS fimage.c fimage.h  flst.c flst.h  flstb_boundary.c flstb_dual.c flstb_dualchain.c flst_bilinear.c flst_boundary.c flstb_quantize.c flst_pixels.c flst_reconstruct.c  fsaddles.c fsignal.c fsignal.h  iio.c iio.h list.c list.h mw_error.c mw_error.h shape.c shape.h)       #Define a variable that contains files of the lib
add_library(FLSTLIB ${FLST_LIB_SRCS}) #Define that those files create the lib FooLib
target_link_libraries(FLSTLIB ${LIBS})

#TODO:Lo mismo hay que poner todos
# ${CMAKE_CURRENT_BINARY_DIR}/iio.o






##setup Config.cmake
#
#
set(PARAMETROS_BASE_DIR "${PROJECT_SOURCE_DIR}")
set(PARAMETROS_CMAKE_DIR "${PROJECT_BINARY_DIR}")
set(FLST_LIBS ${PROJECT_SOURCE_DIR}/libFLSTLIB.a ${LIBS})
set(FLST_INCLUDE ${PROJECT_SOURCE_DIR})
message( "${PROJECT_BINARY_DIR}" )
message( "${PROJECT_SOURCE_DIR}" )
configure_file(flstConfig.cmake.in
   "${PROJECT_BINARY_DIR}/flstConfig.cmake" @ONLY)

