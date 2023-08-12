cmake_minimum_required(VERSION 2.6)
project(point_cloud_feature)
 
find_package(PCL 1.2 REQUIRED)

SET(ENV{VTK_DIR} /usr/include/vtk-7.1)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package (Eigen3 REQUIRED NO_MODULE)
# SET(ENV{EIGEN_INCLUDE_DIRS} /usr/include/eigen3)
# find_package(Eigen REQUIRED)
# find_package(unsopported REQUIRED)
# find_package(eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
 
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
 
add_executable(point_cloud_feature point_cloud_feature.cpp)
 
target_link_libraries (point_cloud_feature ${PCL_LIBRARIES})
 
install(TARGETS point_cloud_feature RUNTIME DESTINATION bin)

