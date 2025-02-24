cmake_minimum_required(VERSION 3.20)
project(root)

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

option(TIME_TESTS, OFF)
option(BUILD_TESTS, OFF)

add_subdirectory(src)
if(TIME_TESTS)
    add_subdirectory(timeMesure)
endif()

if(BUILD_TESTS)
    add_subdirectory(tests)
endif()


# add_executable(root main.cpp)
# target_link_libraries(root PRIVATE LinearSystemSolver)
