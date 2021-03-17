cmake_minimum_required(VERSION 2.6)

project(HOMMGEN)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -std=c++20")

add_subdirectory(GUI)
add_subdirectory(impl)
add_subdirectory(test)