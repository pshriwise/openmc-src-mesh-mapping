cmake_minimum_required(VERSION 3.19)

project(src_mapping)

find_package(OpenMC REQUIRED)

add_executable(omc-src-map-to-mesh src_to_mesh.cpp)
target_link_libraries(omc-src-map-to-mesh OpenMC::libopenmc)
