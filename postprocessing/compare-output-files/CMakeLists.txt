cmake_minimum_required(VERSION 3.8)

project(compare_vtus)

add_executable(${PROJECT_NAME}
   main.cpp
)

# VTK
find_package(VTK)
include_directories(${VTK_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
