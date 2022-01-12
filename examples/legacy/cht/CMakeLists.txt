#find_library(PRECICE_LIB NAMES precice HINTS "${LD_LIBRARY_PATH}")
#find_package(Boost 1.65.1 REQUIRED COMPONENTS log system) #Require same version as preCICE

#string(REPLACE ":" ";" LIBRARY_DIRS $ENV{LD_LIBRARY_PATH})
#find_library(PRECICE_LIB NAMES precice HINTS ${LIBRARY_DIRS})
set(LIBRARY_DIRS "")
if(DEFINED ENV{LD_LIBRARY_PATH})
    message(STATUS "LD_LIBRARY_PATH is set. Use it for hinting")
    string(REPLACE ":" ";" LIBRARY_DIRS $ENV{LD_LIBRARY_PATH} )
endif()
message(STATUS "LIBRARY_DIRS: ${LIBRARY_DIRS}")
find_library(PRECICE_LIB NAMES precice HINTS ${LIBRARY_DIRS})
find_package(Boost 1.65.1 REQUIRED COMPONENTS log system) #Require same version as preCICE

add_subdirectory(monolithic)
add_subdirectory(monolithic-const-dt)
add_subdirectory(partitioned)
#add_subdirectory(iterative-reversed)

