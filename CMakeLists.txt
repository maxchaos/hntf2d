CMAKE_MINIMUM_REQUIRED(VERSION 3.21)
PROJECT(hntf2d VERSION 1.0.0 LANGUAGES CXX)

set(EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 REQUIRED)
find_package(tinyxml2 REQUIRED)

add_library(hntf2d
  SHARED
  cpp/src/map.cpp)

target_include_directories(hntf2d
  PUBLIC
  cpp/include
  ${EIGEN3_INCLUDE_DIRS}
  ${TINYXML2_INCLUDE_DIRS})
