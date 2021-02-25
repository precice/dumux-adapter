cmake_minimum_required(VERSION 3.8)

project(compare_mono_iter_vtu)

add_executable(${PROJECT_NAME}
   main.cpp
)

# VTK
find_package(VTK)
include_directories(${VTK_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})

# Boost
find_package(Boost COMPONENTS filesystem system)
target_include_directories(${PROJECT_NAME} PUBLIC ${BOOST_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
