cmake_minimum_required (VERSION 2.8.11)
project (stereoRectifyMinimal)

add_subdirectory (eigen)

add_executable (stereoRectify stereoRectify.cpp)

target_link_libraries (stereoRectify LINK_PUBLIC eigen)