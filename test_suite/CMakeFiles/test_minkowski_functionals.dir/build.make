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
include test_suite/CMakeFiles/test_minkowski_functionals.dir/depend.make

# Include the progress variables for this target.
include test_suite/CMakeFiles/test_minkowski_functionals.dir/progress.make

# Include the compile flags for this target's objects.
include test_suite/CMakeFiles/test_minkowski_functionals.dir/flags.make

test_suite/CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.o: test_suite/CMakeFiles/test_minkowski_functionals.dir/flags.make
test_suite/CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.o: test_suite/test_minkowski_functionals.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/theorie/spiller/code/ana/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test_suite/CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.o"
	cd /home/theorie/spiller/code/ana/test_suite && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.o -c /home/theorie/spiller/code/ana/test_suite/test_minkowski_functionals.cpp

test_suite/CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.i"
	cd /home/theorie/spiller/code/ana/test_suite && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/theorie/spiller/code/ana/test_suite/test_minkowski_functionals.cpp > CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.i

test_suite/CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.s"
	cd /home/theorie/spiller/code/ana/test_suite && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/theorie/spiller/code/ana/test_suite/test_minkowski_functionals.cpp -o CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.s

# Object files for target test_minkowski_functionals
test_minkowski_functionals_OBJECTS = \
"CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.o"

# External object files for target test_minkowski_functionals
test_minkowski_functionals_EXTERNAL_OBJECTS =

test_suite/bin/test_minkowski_functionals: test_suite/CMakeFiles/test_minkowski_functionals.dir/test_minkowski_functionals.cpp.o
test_suite/bin/test_minkowski_functionals: test_suite/CMakeFiles/test_minkowski_functionals.dir/build.make
test_suite/bin/test_minkowski_functionals: test_suite/CMakeFiles/test_minkowski_functionals.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/theorie/spiller/code/ana/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bin/test_minkowski_functionals"
	cd /home/theorie/spiller/code/ana/test_suite && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_minkowski_functionals.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test_suite/CMakeFiles/test_minkowski_functionals.dir/build: test_suite/bin/test_minkowski_functionals

.PHONY : test_suite/CMakeFiles/test_minkowski_functionals.dir/build

test_suite/CMakeFiles/test_minkowski_functionals.dir/clean:
	cd /home/theorie/spiller/code/ana/test_suite && $(CMAKE_COMMAND) -P CMakeFiles/test_minkowski_functionals.dir/cmake_clean.cmake
.PHONY : test_suite/CMakeFiles/test_minkowski_functionals.dir/clean

test_suite/CMakeFiles/test_minkowski_functionals.dir/depend:
	cd /home/theorie/spiller/code/ana && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/theorie/spiller/code/ana /home/theorie/spiller/code/ana/test_suite /home/theorie/spiller/code/ana /home/theorie/spiller/code/ana/test_suite /home/theorie/spiller/code/ana/test_suite/CMakeFiles/test_minkowski_functionals.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test_suite/CMakeFiles/test_minkowski_functionals.dir/depend

