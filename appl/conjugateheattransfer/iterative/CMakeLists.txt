add_input_file_links()
dune_symlink_to_source_files(FILES "precice-config.xml" "precice-config-serial-implicit.xml")

add_executable(test_freeflow EXCLUDE_FROM_ALL main_freeflow.cc preciceadapter.cc dumuxpreciceindexwrapper.cc)
add_executable(test_solidenergy EXCLUDE_FROM_ALL main_solidenergy.cc preciceadapter.cc dumuxpreciceindexwrapper.cc)
#add_executable(test_preciceadapter EXCLUDE_FROM_ALL preciceadapter.cc)

target_compile_definitions(test_freeflow PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(test_solidenergy PUBLIC "ENABLEMONOLITHIC=0")

#set_property( TARGET libprecice IMPORTED_LOCATION $ENV{PRECICE_ROOT}/build/last/librprecice.so )

find_library(PRECICE_LIB precice)
find_package(Boost 1.65.0 REQUIRED COMPONENTS log system) #Require same version as preCICE
include_directories(${BOOST_INCLUDE_DIR})
add_definitions(-DBOOST_LOG_DYN_LINK)

target_link_libraries(test_freeflow precice ${Boost_LOG_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(test_solidenergy precice ${Boost_LOG_LIBRARY} ${Boost_SYSTEM_LIBRARY})
#target_link_libraries( test_freeflow precice )
#set_property( TARGET libprecice IMPORTED_LOCATION $ENV{PRECICE_ROOT}/build/last/librprecice.so )
