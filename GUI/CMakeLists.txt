cmake_minimum_required(VERSION 2.6)
project(gui)

LIST(APPEND CMAKE_MODULE_PATH "..")

include_directories("..")

add_executable(gui
main.cpp
)

target_link_libraries(gui nana impl Xcursor)

add_custom_command(
        TARGET gui POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_BINARY_DIR}/gui
                ../../bin/gui)
