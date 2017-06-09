cmake_minimum_required(VERSION 2.6)

# DEFAULT BUILD: RELEASE
if (NOT CMAKE_BUILD_TYPE)
   message(STATUS "No build type selected, default to Release")
   set(CMAKE_BUILD_TYPE "Release")
endif()



# IIO (soon to be used)
if( NOT TARGET IIOLIB )
   INCLUDE_DIRECTORIES("3rdparty/iio")
   ADD_SUBDIRECTORY(3rdparty/iio)
ENDIF() 
SET( LIBS ${LIBS} IIOLIB )


 # FLST 
if( NOT TARGET FLSTLIB )
   INCLUDE_DIRECTORIES("3rdparty/flst")
   ADD_SUBDIRECTORY(3rdparty/flst)
ENDIF() 
SET( LIBS ${LIBS} FLSTLIB )



# OpenMP (not used yet)
FIND_PACKAGE( OpenMP )
IF(OPENMP_FOUND)
   MESSAGE("OPENMP FOUND")
   SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}" )
ENDIF()

# IF THE COMPILER IS NOT GCC FAIL!
if(CMAKE_COMPILER_IS_GNUCC)
   SET( EXTRA_COMPILER_FLAGS "${EXTRA_COMPILER_FLAGS} -Wall" )
   # compiler flags for 3rdparty/simpois
   set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else()
   MESSAGE("This project needs c99")
endif()

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()

# UNCOMMENT TO ENABLE DEBUG OUTPUT
#SET( EXTRA_COMPILER_FLAGS "${EXTRA_COMPILER_FLAGS} -DDBG_OUTPUT" )

# SET THE FLAG FOR THE IPOL DEMO 
if (IPOL_DEMO)
   message(STATUS "Generating the exacutable for the IPOL demo")
   SET( EXTRA_COMPILER_FLAGS "${EXTRA_COMPILER_FLAGS} -DIPOL_DEMO" )
endif()


add_executable (FLST_test   
                point.cpp
                distance_transform.cpp                  
                gaussian_weights.cpp  
                bilateral_filter.cpp      
                main_flst.cpp                    
                io_utility.cpp 
                flst_utility.cpp             
                mask_iterator.cpp           
                patchsize.cpp
                gradient.cpp
                mask.cpp
                image.hpp
                image.h       
                point.h
                distance_transform.h                
                gaussian_weights.h  
                bilateral_filter.h      
                io_utility.h    
                flst_utility.h          
                mask_iterator.h           
                patchsize.h
                gradient.h
                mask.h
                3rdparty/simpois/simpois.c
               )
target_link_libraries(FLST_test ${LIBS})
SET_TARGET_PROPERTIES(FLST_test PROPERTIES
      COMPILE_FLAGS "${EXTRA_COMPILER_FLAGS} ${OpenMP_CXX_FLAGS}")