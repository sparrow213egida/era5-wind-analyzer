cmake_minimum_required(VERSION 4.0)
project(TestTask)

set(CMAKE_CXX_STANDARD 20)

find_package(GDAL REQUIRED)

add_executable(TestTask main.cpp)

target_link_libraries(TestTask PRIVATE GDAL::GDAL)
