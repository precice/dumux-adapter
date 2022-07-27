# Monolithic test cases

# 2D cases
add_subdirectory(td-ff-pm-2d-chaabani) # Chaabani et al.

# 3D cases
if ( dune-subgrid_FOUND )
    add_subdirectory(flow-over-step-3d) # Chaabani et al.
endif()