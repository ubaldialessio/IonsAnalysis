# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /cvmfs/sft.cern.ch/lcg/releases/CMake/3.18.4-2ffec/x86_64-centos7-gcc9-opt/bin/cmake

# The command to remove a file.
RM = /cvmfs/sft.cern.ch/lcg/releases/CMake/3.18.4-2ffec/x86_64-centos7-gcc9-opt/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /storage/gpfs_ams/ams/users/aubaldi/checkVar

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /storage/gpfs_ams/ams/users/aubaldi/checkVar/build

# Include any dependencies generated for this target.
include CMakeFiles/checkVar.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/checkVar.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/checkVar.dir/flags.make

CMakeFiles/checkVar.dir/src/checkVar.cpp.o: CMakeFiles/checkVar.dir/flags.make
CMakeFiles/checkVar.dir/src/checkVar.cpp.o: ../src/checkVar.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/storage/gpfs_ams/ams/users/aubaldi/checkVar/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/checkVar.dir/src/checkVar.cpp.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/9.3.0-6991b/x86_64-centos7/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/checkVar.dir/src/checkVar.cpp.o -c /storage/gpfs_ams/ams/users/aubaldi/checkVar/src/checkVar.cpp

CMakeFiles/checkVar.dir/src/checkVar.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/checkVar.dir/src/checkVar.cpp.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/9.3.0-6991b/x86_64-centos7/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /storage/gpfs_ams/ams/users/aubaldi/checkVar/src/checkVar.cpp > CMakeFiles/checkVar.dir/src/checkVar.cpp.i

CMakeFiles/checkVar.dir/src/checkVar.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/checkVar.dir/src/checkVar.cpp.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/9.3.0-6991b/x86_64-centos7/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /storage/gpfs_ams/ams/users/aubaldi/checkVar/src/checkVar.cpp -o CMakeFiles/checkVar.dir/src/checkVar.cpp.s

# Object files for target checkVar
checkVar_OBJECTS = \
"CMakeFiles/checkVar.dir/src/checkVar.cpp.o"

# External object files for target checkVar
checkVar_EXTERNAL_OBJECTS =

../bin/checkVar: CMakeFiles/checkVar.dir/src/checkVar.cpp.o
../bin/checkVar: CMakeFiles/checkVar.dir/build.make
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/naia/v1.0.2/lib/libNAIAChain.so
../bin/checkVar: /storage/gpfs_ams/ams/users/aubaldi/nsl.install/lib/libNSLSelections.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/naia/v1.0.2/lib/libNAIAContainers.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/naia/v1.0.2/lib/libNAIAUtility.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libCore.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libImt.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libRIO.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libNet.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libHist.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libGraf.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libGraf3d.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libGpad.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libROOTDataFrame.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libTree.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libTreePlayer.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libRint.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libPostscript.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libMatrix.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libPhysics.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libMathCore.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libThread.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libMultiProc.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libROOTVecOps.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/naia/v1.0.2/lib64/libspdlog.a
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/naia/v1.0.2/lib64/libfmt.a
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libPhysics.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libPostscript.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libROOTDataFrame.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libROOTVecOps.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libRint.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libTreePlayer.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libGraf3d.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libGpad.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libGraf.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libHist.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libMatrix.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libTree.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libMathCore.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libImt.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libMultiProc.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libNet.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libRIO.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libThread.so
../bin/checkVar: /cvmfs/ams.cern.ch/Offline/amsitaly/public/install/x86_64-centos7-gcc9.3/root-6.26.02/lib/libCore.so
../bin/checkVar: CMakeFiles/checkVar.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/storage/gpfs_ams/ams/users/aubaldi/checkVar/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/checkVar"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/checkVar.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/checkVar.dir/build: ../bin/checkVar

.PHONY : CMakeFiles/checkVar.dir/build

CMakeFiles/checkVar.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/checkVar.dir/cmake_clean.cmake
.PHONY : CMakeFiles/checkVar.dir/clean

CMakeFiles/checkVar.dir/depend:
	cd /storage/gpfs_ams/ams/users/aubaldi/checkVar/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /storage/gpfs_ams/ams/users/aubaldi/checkVar /storage/gpfs_ams/ams/users/aubaldi/checkVar /storage/gpfs_ams/ams/users/aubaldi/checkVar/build /storage/gpfs_ams/ams/users/aubaldi/checkVar/build /storage/gpfs_ams/ams/users/aubaldi/checkVar/build/CMakeFiles/checkVar.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/checkVar.dir/depend
