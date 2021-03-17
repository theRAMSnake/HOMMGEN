cmake_minimum_required(VERSION 2.6)
project(impl)

LIST(APPEND CMAKE_MODULE_PATH "..")

add_library(
impl 
generator.cpp
rng.cpp
objectlib.cpp
monsterlib.cpp
)

find_package(Boost COMPONENTS system REQUIRED)

#target_link_libraries(impl )
