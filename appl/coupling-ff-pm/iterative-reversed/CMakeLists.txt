# executables for the iterative case
string(REPLACE ":" ";" LIBRARY_DIRS $ENV{LD_LIBRARY_PATH})
find_library(PRECICE_LIB NAMES precice HINTS ${LIBRARY_DIRS})
find_package(Boost 1.65.1 REQUIRED COMPONENTS log system) #Require same version as preCICE

#dune_add_test(NAME test_ff_reversed
#              SOURCES main_ff-reversed.cc ../precice/preciceadapter.cc ../precice/dumuxpreciceindexwrapper.cc
#              COMPILE_DEFINITIONS ENABLEMONOLITHIC=0)

#dune_add_test(NAME test_pm_reversed
#              SOURCES main_pm-reversed.cc ../precice/preciceadapter.cc ../precice/dumuxpreciceindexwrapper.cc
#              COMPILE_DEFINITIONS ENABLEMONOLITHIC=0)

add_executable(test_ff_reversed EXCLUDE_FROM_ALL main_ff-reversed.cc ../../precice-adapter/src/preciceadapter.cc ../../precice-adapter/src/dumuxpreciceindexwrapper.cc)
add_executable(test_pm_reversed EXCLUDE_FROM_ALL main_pm-reversed.cc ../../precice-adapter/src/preciceadapter.cc ../../precice-adapter/src/dumuxpreciceindexwrapper.cc)

target_compile_definitions(test_ff_reversed PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(test_pm_reversed PUBLIC "ENABLEMONOLITHIC=0")

target_compile_definitions(test_ff_reversed PRIVATE BOOST_ALL_DYN_LINK)
target_link_libraries(test_ff_reversed PRIVATE ${Boost_LIBRARIES} ${PRECICE_LIB})
target_include_directories(test_ff_reversed PRIVATE ${Boost_INCLUDE_DIRS})

target_compile_definitions(test_pm_reversed PRIVATE BOOST_ALL_DYN_LINK)
target_link_libraries(test_pm_reversed PRIVATE ${Boost_LIBRARIES} ${PRECICE_LIB})
target_include_directories(test_pm_reversed PRIVATE ${Boost_INCLUDE_DIRS})
# add a symlink for each input file
add_input_file_links()
