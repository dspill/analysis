# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /usr/data/spiller/local/bin/cmake

# The command to remove a file.
RM = /usr/data/spiller/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/theorie/spiller/code/ana

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/theorie/spiller/code/ana

# Include any dependencies generated for this target.
include src/CMakeFiles/get_last_frame.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/get_last_frame.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/get_last_frame.dir/flags.make

src/CMakeFiles/get_last_frame.dir/get_last_frame.cpp.o: src/CMakeFiles/get_last_frame.dir/flags.make
src/CMakeFiles/get_last_frame.dir/get_last_frame.cpp.o: src/get_last_frame.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/theorie/spiller/code/ana/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/get_last_frame.dir/get_last_frame.cpp.o"
	cd /home/theorie/spiller/code/ana/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/get_last_frame.dir/get_last_frame.cpp.o -c /home/theorie/spiller/code/ana/src/get_last_frame.cpp

src/CMakeFiles/get_last_frame.dir/get_last_frame.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/get_last_frame.dir/get_last_frame.cpp.i"
	cd /home/theorie/spiller/code/ana/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/theorie/spiller/code/ana/src/get_last_frame.cpp > CMakeFiles/get_last_frame.dir/get_last_frame.cpp.i

src/CMakeFiles/get_last_frame.dir/get_last_frame.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/get_last_frame.dir/get_last_frame.cpp.s"
	cd /home/theorie/spiller/code/ana/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/theorie/spiller/code/ana/src/get_last_frame.cpp -o CMakeFiles/get_last_frame.dir/get_last_frame.cpp.s

# Object files for target get_last_frame
get_last_frame_OBJECTS = \
"CMakeFiles/get_last_frame.dir/get_last_frame.cpp.o"

# External object files for target get_last_frame
get_last_frame_EXTERNAL_OBJECTS =

bin/get_last_frame: src/CMakeFiles/get_last_frame.dir/get_last_frame.cpp.o
bin/get_last_frame: src/CMakeFiles/get_last_frame.dir/build.make
bin/get_last_frame: src/CMakeFiles/get_last_frame.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/theorie/spiller/code/ana/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/get_last_frame"
	cd /home/theorie/spiller/code/ana/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/get_last_frame.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/get_last_frame.dir/build: bin/get_last_frame

.PHONY : src/CMakeFiles/get_last_frame.dir/build

src/CMakeFiles/get_last_frame.dir/clean:
	cd /home/theorie/spiller/code/ana/src && $(CMAKE_COMMAND) -P CMakeFiles/get_last_frame.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/get_last_frame.dir/clean

src/CMakeFiles/get_last_frame.dir/depend:
	cd /home/theorie/spiller/code/ana && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/theorie/spiller/code/ana /home/theorie/spiller/code/ana/src /home/theorie/spiller/code/ana /home/theorie/spiller/code/ana/src /home/theorie/spiller/code/ana/src/CMakeFiles/get_last_frame.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/get_last_frame.dir/depend

