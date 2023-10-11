cmake_minimum_required(VERSION 3.19)

project(src_mapping)

find_package(OpenMC REQUIRED)
find_package(OpenMP)

add_executable(omc-src-map-to-mesh src_to_mesh.cpp)
target_link_libraries(omc-src-map-to-mesh OpenMC::libopenmc)

if (OpenMP_FOUND)
    target_compile_definitions(omc-src-map-to-mesh PRIVATE -DSRC_MAP_OPENMP)
    target_link_libraries(omc-src-map-to-mesh OpenMP::OpenMP_CXX)
endif()