add_input_file_links()
dune_symlink_to_source_files(FILES "precice-config.xml" "precice-config-serial-implicit.xml")

add_executable(test_freeflow EXCLUDE_FROM_ALL main_freeflow.cc ../common/preciceadapter.cc ../common/dumuxpreciceindexwrapper.cc ../common/helperfunctions.cc)
add_executable(test_solidenergy EXCLUDE_FROM_ALL main_solidenergy.cc ../common/preciceadapter.cc ../common/dumuxpreciceindexwrapper.cc ../common/helperfunctions.cc)
#add_executable(test_preciceadapter EXCLUDE_FROM_ALL preciceadapter.cc)

target_compile_definitions(test_freeflow PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(test_solidenergy PUBLIC "ENABLEMONOLITHIC=0")

#set_property( TARGET libprecice IMPORTED_LOCATION $ENV{PRECICE_ROOT}/build/last/librprecice.so )

string(REPLACE ":" ";" LIBRARY_DIRS $ENV{LD_LIBRARY_PATH})
find_library(PRECICE_LIB NAMES precice HINTS "${LIBRARY_DIRS}")
find_package(Boost 1.65.1 REQUIRED COMPONENTS log system) #Require same version as preCICE


target_compile_definitions(test_freeflow PRIVATE BOOST_ALL_DYN_LINK)
target_link_libraries(test_freeflow PRIVATE ${Boost_LIBRARIES} ${PRECICE_LIB})
target_include_directories(test_freeflow PRIVATE ${Boost_INCLUDE_DIRS})

target_compile_definitions(test_solidenergy PRIVATE BOOST_ALL_DYN_LINK)
target_link_libraries(test_solidenergy PRIVATE ${Boost_LIBRARIES} ${PRECICE_LIB})
target_include_directories(test_solidenergy PRIVATE ${Boost_INCLUDE_DIRS})

#include_directories(${BOOST_INCLUDE_DIR})
#add_definitions(-DBOOST_LOG_DYN_LINK)

#target_link_libraries(test_freeflow precice ${Boost_LOG_LIBRARY} ${Boost_SYSTEM_LIBRARY})
#target_link_libraries(test_solidenergy precice ${Boost_LOG_LIBRARY} ${Boost_SYSTEM_LIBRARY})
#target_link_libraries( test_freeflow precice )
#set_property( TARGET libprecice IMPORTED_LOCATION $ENV{PRECICE_ROOT}/build/last/librprecice.so )
