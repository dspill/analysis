add_library(Real3D     Real3D.hpp)
add_library(Vector     vector.hpp)
add_library(Couf       couf.hpp)

add_library(Trajectory trajectory.hpp)
add_library(Frame      frame.hpp)
add_library(Lb_frame   lb_frame.hpp)
add_library(Molecule   molecule.hpp)
add_library(Timeseries timeseries.hpp)

# Make sure the compiler can find include files for our libraries
# when other libraries or executables link to them
target_include_directories(Real3D     PUBLIC ${CMAKE_CURRWNT_SOURCE_DIR})
target_include_directories(Vector     PUBLIC ${CMAKE_CURRWNT_SOURCE_DIR})
target_include_directories(Couf       PUBLIC ${CMAKE_CURRWNT_SOURCE_DIR})

target_include_directories(Trajectory PUBLIC ${CMAKE_CURRWNT_SOURCE_DIR})
target_include_directories(Frame      PUBLIC ${CMAKE_CURRWNT_SOURCE_DIR})
target_include_directories(Lb_frame   PUBLIC ${CMAKE_CURRWNT_SOURCE_DIR})
target_include_directories(Molecule   PUBLIC ${CMAKE_CURRWNT_SOURCE_DIR})
target_include_directories(Timeseries PUBLIC ${CMAKE_CURRWNT_SOURCE_DIR})
