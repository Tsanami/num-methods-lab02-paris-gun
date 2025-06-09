cmake_minimum_required(VERSION 3.10)
project(ParishGunModel CXX)
set(CMAKE_CXX_STANDARD 17)
add_library(geodesy geodesy.cpp)
add_library(atmosphere atmosphere.cpp)
add_executable(run_test main.cpp)
target_link_libraries(run_test geodesy atmosphere)

