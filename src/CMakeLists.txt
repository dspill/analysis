set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

set(programs 
    analyze
    analyze_lb
    analyze_timeseries 
    analyze_trajectory 
    average_timeseries 
    brute_force_dynamic_structure_factor 
    check
    correlation_function 
    dynamic_structure_factor 
    dynamic_structure_factor_2d
    dynamic_structure_factor_from_lattice
    generate_timeseries_ete_over_rg 
    get_frames
    get_last_frame 
    latticeize_configuration 
    manipulate_frame 
    minkowski_functionals
    r_e_over_r_g
    slice_configuration
    traj_info
    ts_correlation_function 
    concatenate_trajectories 
    )

foreach(prog ${programs})
    add_executable(${prog} ${prog}.cpp)
    target_link_libraries(${prog} fftw3)
endforeach()
