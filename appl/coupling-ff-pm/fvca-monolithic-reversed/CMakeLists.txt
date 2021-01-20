
add_input_file_links()

add_executable(fvca-monolithic-reversed EXCLUDE_FROM_ALL fvca-monolithic.cc)

dune_add_test(NAME test_fvca_monolithic_reversed
              TARGET fvca-monolithic-reversed
              LABELS ff-pm stokes monolithic darcy
              TIMEOUT 10
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/monolithic-new/stokes-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/stokes-darcy-monolithic-coupling-test_stokes-00001.vtu
                      ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/monolithic-new/darcy-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/stokes-darcy-monolithic-coupling-test_darcy-00001.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/fvca-monolithic-reversed ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/monolithic-new/params.input"
              --relative 1e-5
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)