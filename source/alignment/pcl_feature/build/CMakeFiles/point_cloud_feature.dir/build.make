# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature/build

# Include any dependencies generated for this target.
include CMakeFiles/point_cloud_feature.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/point_cloud_feature.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/point_cloud_feature.dir/flags.make

CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.o: CMakeFiles/point_cloud_feature.dir/flags.make
CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.o: ../point_cloud_feature.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.o -c /home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature/point_cloud_feature.cpp

CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature/point_cloud_feature.cpp > CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.i

CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature/point_cloud_feature.cpp -o CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.s

# Object files for target point_cloud_feature
point_cloud_feature_OBJECTS = \
"CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.o"

# External object files for target point_cloud_feature
point_cloud_feature_EXTERNAL_OBJECTS =

point_cloud_feature: CMakeFiles/point_cloud_feature.dir/point_cloud_feature.cpp.o
point_cloud_feature: CMakeFiles/point_cloud_feature.dir/build.make
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_apps.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_outofcore.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_people.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libboost_system.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libboost_regex.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libqhull.so
point_cloud_feature: /usr/lib/libOpenNI.so
point_cloud_feature: /usr/lib/libOpenNI2.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkChartsCore-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkInfovisCore-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libfreetype.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libz.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libjpeg.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpng.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libtiff.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libexpat.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkIOGeometry-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkIOLegacy-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkIOPLY-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkRenderingLOD-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkViewsContext2D-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkViewsCore-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkRenderingContextOpenGL2-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkRenderingOpenGL2-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libflann_cpp.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_surface.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_keypoints.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_tracking.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_recognition.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_registration.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_stereo.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_segmentation.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_features.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_filters.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_sample_consensus.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_ml.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_visualization.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_search.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_kdtree.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_io.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_octree.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libpcl_common.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkInteractionWidgets-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkFiltersModeling-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkInteractionStyle-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkFiltersExtraction-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkFiltersStatistics-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkImagingFourier-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkalglib-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkFiltersHybrid-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkImagingGeneral-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkImagingSources-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkImagingHybrid-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkRenderingAnnotation-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkImagingColor-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkRenderingVolume-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkIOXML-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkIOXMLParser-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkIOCore-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkRenderingContext2D-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkRenderingFreeType-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libfreetype.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkImagingCore-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkRenderingCore-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkCommonColor-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkFiltersGeometry-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkFiltersSources-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkFiltersGeneral-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkCommonComputationalGeometry-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkFiltersCore-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkIOImage-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkCommonExecutionModel-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkCommonDataModel-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkCommonTransforms-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkCommonMisc-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkCommonMath-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkCommonSystem-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkCommonCore-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtksys-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkDICOMParser-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libvtkmetaio-7.1.so.7.1p.1
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libz.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libGLEW.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libSM.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libICE.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libX11.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libXext.so
point_cloud_feature: /usr/lib/x86_64-linux-gnu/libXt.so
point_cloud_feature: CMakeFiles/point_cloud_feature.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable point_cloud_feature"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/point_cloud_feature.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/point_cloud_feature.dir/build: point_cloud_feature

.PHONY : CMakeFiles/point_cloud_feature.dir/build

CMakeFiles/point_cloud_feature.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/point_cloud_feature.dir/cmake_clean.cmake
.PHONY : CMakeFiles/point_cloud_feature.dir/clean

CMakeFiles/point_cloud_feature.dir/depend:
	cd /home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature /home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature /home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature/build /home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature/build /home/data/sleeping_work/Structure-retrieval/em_retrieval/github_repo/source/alignment/pcl_feature/build/CMakeFiles/point_cloud_feature.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/point_cloud_feature.dir/depend

