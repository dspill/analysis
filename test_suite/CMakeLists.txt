set (CMAKE_BUILD_TYPE Debug)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")

set(test_programs frame minkowski molecule mctimeseries timeseries trajectory vector)

foreach(tprog ${test_programs})
    add_executable(${tprog} ${tprog}.cpp)
    target_link_libraries(${tprog} fftw3)
endforeach()

if(NOT TARGET test_suite)
  add_custom_target(test_suite)
  add_dependencies(test_suite ${test_programs})
endif()
