cmake_minimum_required(VERSION 3.19)
project(FO2)

set(CMAKE_CXX_STANDARD 20)

#FILE(GLOB_RECURSE files *.cpp *.hpp)

add_definitions(-DPROJECT_DIR_RAW="${PROJECT_SOURCE_DIR}")

#include_directories(/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3)
#include_directories(/usr/include/eigen3)
include_directories(/home/xidad/libs/eigen-3.4.0/)
#add_compile_options("-O3")


add_executable(FO2 main.cpp utils.hpp)
