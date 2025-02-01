cmake_minimum_required(VERSION 3.20)
project(root)

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")


add_subdirectory(src)

add_executable(root main.cpp)
target_link_libraries(root PRIVATE LinearSystemSolver)

add_subdirectory(tests)