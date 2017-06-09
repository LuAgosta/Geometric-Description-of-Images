FILE(REMOVE_RECURSE
  "iioConfig.cmake"
  "CMakeCache.txt"
  "CMakeFiles"
  "Makefile"
  "cmake_install.cmake"
  "iio"
  "iion"
  "iio.o"
  "CMakeFiles/IIOLIB.dir/iio.c.o"
  "libIIOLIB.pdb"
  "libIIOLIB.a"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang C)
  INCLUDE(CMakeFiles/IIOLIB.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
