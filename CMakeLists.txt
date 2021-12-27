cmake_minimum_required(VERSION 3.1)
project(FindInnerBoundary)


include_directories(include)

link_directories("lib")


add_executable(FindInnerBoundary 
code/main.cpp
code/HoleFiller.h
code/HoleFiller.cpp
code/MeshDoctor.h
code/MeshDoctor.cpp
code/FindInnerBoundary.h
code/FindInnerBoundary.cpp
)

target_link_libraries(FindInnerBoundary debug OpenMeshCored.lib)
target_link_libraries(FindInnerBoundary optimized OpenMeshCore.lib)