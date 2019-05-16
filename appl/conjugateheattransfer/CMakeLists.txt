#find_library(PRECICE_LIB NAMES precice HINTS "${LD_LIBRARY_PATH}")
#find_package(Boost 1.65.1 REQUIRED COMPONENTS log system) #Require same version as preCICE

string(REPLACE ":" ";" LIBRARY_DIRS $ENV{LD_LIBRARY_PATH})
find_library(PRECICE_LIB NAMES precice HINTS ${LIBRARY_DIRS})
find_package(Boost 1.65.1 REQUIRED COMPONENTS log system) #Require same version as preCICE

add_subdirectory(monolithic)
add_subdirectory(monolithic-const-dt)
add_subdirectory(iterative)
#add_subdirectory(iterative-reversed)

