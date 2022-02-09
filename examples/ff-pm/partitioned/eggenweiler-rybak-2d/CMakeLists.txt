add_input_file_links()
#dune_symlink_to_source_files(FILES "convergencetest.py")

#add_executable(test_md_boundary_darcy1ptpfa_freeflow1p_convtest EXCLUDE_FROM_ALL main.cc)
#target_compile_definitions(test_md_boundary_darcy1ptpfa_freeflow1p_convtest PUBLIC DARCYTYPETAG=DarcyOnePCC)

add_executable(test_eggenweiler_rybak_2d_darcy EXCLUDE_FROM_ALL main_darcy.cc)
target_link_libraries(test_eggenweiler_rybak_2d_darcy PRIVATE dumux-precice)
#target_compile_definitions(test_md_boundary_darcy1pbox_freeflow1p_convtest PUBLIC DARCYTYPETAG=DarcyOnePBox)

add_executable(test_eggenweiler_rybak_2d_stokes EXCLUDE_FROM_ALL main_stokes.cc)
target_link_libraries(test_eggenweiler_rybak_2d_stokes PRIVATE dumux-precice)
#target_compile_definitions(test_md_boundary_darcy1pbox_freeflow1p_convtest PUBLIC DARCYTYPETAG=DarcyOnePBox)


dune_add_test(NAME test_md_boundary_darcy1p_stokes1p_newic
              TARGET test_eggenweiler_rybak_2d_darcy test_eggenweiler_rybak_2d_stokes
              LABELS multidomain multidomain_boundary stokesdarcy eggenweilerrybak
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS -e test_md_boundary_darcy1pbox_freeflow1p_convtest
                       -i params.input
                       -c EggenweilerRybak
                       -Problem.SlipCondition ER
                       -Problem.NewIc true
                       -FreeFlow.EnableUnsymmetrizedVelocityGradientForIC true
                       -FreeFlow.EnableUnsymmetrizedVelocityGradient true
                       -Darcy.SpatialParams.Permeability "1.0 1.0 0.0"
                       -Darcy.InterfaceParams.EpsInterface 1.0
                       -Darcy.InterfaceParams.N_s_bl 0.5
                       -Darcy.InterfaceParams.N_1_bl -3.0
                       -Darcy.InterfaceParams.M_bl "1.0 1.0"
                       -Darcy.SpatialParams.Porosity 1.0
                       -Component.LiquidKinematicViscosity 2.0)


dune_add_test(NAME dummy_test_md_boundary_darcy1p_stokes1p_newic
                       TARGET test_eggenweiler_rybak_2d_stokes
                       LABELS multidomain multidomain_boundary stokesdarcy eggenweilerrybak dummy
                       TIMEOUT 10
                       CMAKE_GUARD HAVE_UMFPACK
                       COMMAND echo
                       CMD_ARGS "dummy_test_md_boundary_darcy1p_stokes1p_newic")