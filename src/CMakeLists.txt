set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")

set(programs 
    traj_info
    analyze_trajectory 
    analyze_timeseries 
    analyze_lb
    r_e_over_r_g
    correlation_function 
    ts_correlation_function 
    generate_timeseries_ete_over_rg 
    minkowski_functionals
    dynamic_structure_factor 
    dynamic_structure_factor_from_lattice
    brute_force_dynamic_structure_factor 
    slice_configuration
    latticeize_configuration 
    manipulate_frame 
    get_last_frame 
    get_frames
    average_timeseries 
    )

foreach(prog ${programs})
    add_executable(${prog} ${prog}.cpp)
endforeach()

target_link_libraries(dynamic_structure_factor_from_lattice fftw3)
target_link_libraries(dynamic_structure_factor fftw3)
