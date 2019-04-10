cmake_minimum_required(VERSION 3.11)
project("ana" C CXX)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wshadow -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
#set (CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fno-omit-frame-pointer -fsanitize=address")
#set (CMAKE_LINKER_FLAGS_DEBUG "-fno-omit-frame-pointer -fsanitize=address")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

include_directories(lib)
add_subdirectory(src)
add_subdirectory(test_suite)

enable_testing()
set(CTEST_OUTPUT_ON_FAILURE 1)
set(testbin test_suite/bin)
file(MAKE_DIRECTORY ${testbin})
set(test_programs frame minkowski molecule timeseries trajectory vector)

foreach(tprog ${test_programs})
    add_test(compile_${tprog} make ${tprog})
    add_test(${tprog} ${testbin}/${tprog})
endforeach()

