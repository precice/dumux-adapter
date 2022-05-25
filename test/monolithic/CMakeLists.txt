# 2D cases
add_subdirectory(flow-over-square-2d)
if ( dune-subgrid_FOUND )
    add_subdirectory(flow-over-step-2d) # Chaabani et al.
endif()
# 3D cases
add_subdirectory(ff-pm-3d)