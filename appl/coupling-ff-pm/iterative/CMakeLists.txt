# executables for the iterative case
set(LIBRARY_DIRS "")
if(DEFINED ENV{LD_LIBRARY_PATH})
    message(STATUS "LD_LIBRARY_PATH is set. Use it for hinting")
    string(REPLACE ":" ";" LIBRARY_DIRS $ENV{LD_LIBRARY_PATH} )
endif()
message(STATUS "LIBRARY_DIRS: ${LIBRARY_DIRS}")
find_library(PRECICE_LIB NAMES precice HINTS ${LIBRARY_DIRS})

find_package(Boost 1.65.1 REQUIRED COMPONENTS log system) #Require same version as preCICE

#dune_add_test(NAME test_ff
#              SOURCES main_ff.cc ../precice/preciceadapter.cc ../precice/dumuxpreciceindexwrapper.cc
#              COMPILE_DEFINITIONS ENABLEMONOLITHIC=0)

#dune_add_test(NAME test_pm
#              SOURCES main_pm.cc ../precice/preciceadapter.cc ../precice/dumuxpreciceindexwrapper.cc
#              COMPILE_DEFINITIONS ENABLEMONOLITHIC=0)

add_executable(test_ff EXCLUDE_FROM_ALL main_ff.cc ../../precice-adapter/src/preciceadapter.cc ../../precice-adapter/src/dumuxpreciceindexwrapper.cc)
add_executable(test_pm EXCLUDE_FROM_ALL main_pm.cc ../../precice-adapter/src/preciceadapter.cc ../../precice-adapter/src/dumuxpreciceindexwrapper.cc)

target_compile_definitions(test_ff PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(test_pm PUBLIC "ENABLEMONOLITHIC=0")

target_compile_definitions(test_ff PRIVATE BOOST_ALL_DYN_LINK)
target_link_libraries(test_ff PRIVATE ${Boost_LIBRARIES} ${PRECICE_LIB})
target_include_directories(test_ff PRIVATE ${Boost_INCLUDE_DIRS})

target_compile_definitions(test_pm PRIVATE BOOST_ALL_DYN_LINK)
target_link_libraries(test_pm PRIVATE ${Boost_LIBRARIES} ${PRECICE_LIB})
target_include_directories(test_pm PRIVATE ${Boost_INCLUDE_DIRS})
# add a symlink for each input file
add_input_file_links()
