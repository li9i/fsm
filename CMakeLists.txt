cmake_minimum_required(VERSION 2.8.3)
project(smsm_tests)

#https://stackoverflow.com/a/38511897
SET(CMAKE_CXX_COMPILER /usr/bin/g++)
add_compile_options(-std=c++11 -Ofast)

find_package(PkgConfig)
pkg_check_modules(csm REQUIRED csm)

link_directories(${csm_LIBRARY_DIRS})

find_package(CGAL REQUIRED)

INCLUDE_DIRECTORIES(include
  ${csm_INCLUDE_DIRS}
  /usr/include/eigen3/)

# smsm #########################################################################
add_library(sm_lib src/sm.cpp)

target_link_libraries(sm_lib fftw3 CGAL::CGAL_Core)

# The executable ###############################################################
add_executable(sm_node
  src/sm_node.cpp)

target_link_libraries(sm_node sm_lib ${csm_LIBRARIES})
