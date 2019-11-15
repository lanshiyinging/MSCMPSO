cmake_minimum_required(VERSION 3.8)
project(MSCMPSO)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp function.h function.cpp parameter.h Macmpso.cpp Macmpso.h)
add_executable(MSCMPSO ${SOURCE_FILES})