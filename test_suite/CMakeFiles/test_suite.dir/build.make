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

# Utility rule file for test_suite.

# Include the progress variables for this target.
include test_suite/CMakeFiles/test_suite.dir/progress.make

test_suite: test_suite/CMakeFiles/test_suite.dir/build.make

.PHONY : test_suite

# Rule to build all files generated by this target.
test_suite/CMakeFiles/test_suite.dir/build: test_suite

.PHONY : test_suite/CMakeFiles/test_suite.dir/build

test_suite/CMakeFiles/test_suite.dir/clean:
	cd /home/theorie/spiller/code/ana/test_suite && $(CMAKE_COMMAND) -P CMakeFiles/test_suite.dir/cmake_clean.cmake
.PHONY : test_suite/CMakeFiles/test_suite.dir/clean

test_suite/CMakeFiles/test_suite.dir/depend:
	cd /home/theorie/spiller/code/ana && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/theorie/spiller/code/ana /home/theorie/spiller/code/ana/test_suite /home/theorie/spiller/code/ana /home/theorie/spiller/code/ana/test_suite /home/theorie/spiller/code/ana/test_suite/CMakeFiles/test_suite.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test_suite/CMakeFiles/test_suite.dir/depend

