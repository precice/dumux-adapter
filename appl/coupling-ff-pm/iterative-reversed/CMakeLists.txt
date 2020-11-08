# executables for the iterative case
string(REPLACE ":" ";" LIBRARY_DIRS $ENV{LD_LIBRARY_PATH})
find_library(PRECICE_LIB NAMES precice HINTS ${LIBRARY_DIRS})
find_package(Boost 1.65.1 REQUIRED COMPONENTS log system) #Require same version as preCICE

#dune_add_test(NAME test_ff_reversed
#              SOURCES main_ff-reversed.cc ../precice/preciceadapter.cc ../precice/dumuxpreciceindexwrapper.cc
#              COMPILE_DEFINITIONS ENABLEMONOLITHIC=0)

#dune_add_test(NAME test_pm_reversed
#              SOURCES main_pm-reversed.cc ../precice/preciceadapter.cc ../precice/dumuxpreciceindexwrapper.cc
#              COMPILE_DEFINITIONS ENABLEMONOLITHIC=0)

add_executable(test_ff_reversed EXCLUDE_FROM_ALL main_ff-reversed.cc ../../precice-adapter/src/preciceadapter.cc ../../precice-adapter/src/dumuxpreciceindexwrapper.cc)
add_executable(test_pm_reversed EXCLUDE_FROM_ALL main_pm-reversed.cc ../../precice-adapter/src/preciceadapter.cc ../../precice-adapter/src/dumuxpreciceindexwrapper.cc)

target_compile_definitions(test_ff_reversed PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(test_pm_reversed PUBLIC "ENABLEMONOLITHIC=0")

target_compile_definitions(test_ff_reversed PRIVATE BOOST_ALL_DYN_LINK)
target_link_libraries(test_ff_reversed PRIVATE ${Boost_LIBRARIES} ${PRECICE_LIB})
target_include_directories(test_ff_reversed PRIVATE ${Boost_INCLUDE_DIRS})

target_compile_definitions(test_pm_reversed PRIVATE BOOST_ALL_DYN_LINK)
target_link_libraries(test_pm_reversed PRIVATE ${Boost_LIBRARIES} ${PRECICE_LIB})
target_include_directories(test_pm_reversed PRIVATE ${Boost_INCLUDE_DIRS})
# add a symlink for each input file
add_input_file_links()

macro(add_precice_file_links)
  FILE(GLOB precice_files *.xml)
  foreach(VAR ${input_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()
add_precice_file_links()



dune_symlink_to_source_files(FILES run-iterative-simulation-test.sh)

#dumux_add_test(NAME test_ff_pm_coupling_serial_implicit_stokes_first
#              TARGET test_ff_reversed test_pm_reversed
#              LABELS freeflow stokes precice darcy
#              CMAKE_GUARD HAVE_UMFPACK
#              COMMAND run-iterative-simulation-test.sh
#              CMD_ARGS ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/serial-implicit/stokes-first/params.input ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/serial-implicit/stokes-first/precice-config.xml)

#dumux_add_test(NAME test_ff_pm_coupling_serial_implicit_darcy_first
#              TARGET test_ff_reversed test_pm_reversed
#              LABELS freeflow stokes precice darcy
#              CMAKE_GUARD HAVE_UMFPACK
#              COMMAND run-iterative-simulation-test.sh
#              CMD_ARGS ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/serial-implicit/darcy-first/params.input ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/serial-implicit/darcy-first/precice-config.xml)

#dumux_add_test(NAME test_ff_pm_coupling_parallel_implicit
#              TARGET test_ff_reversed test_pm_reversed
#              LABELS freeflow stokes precice darcy
#              CMAKE_GUARD HAVE_UMFPACK
#              COMMAND run-iterative-simulation-test.sh
#              CMD_ARGS ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/params.input ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/precice-config.xml)              


dumux_add_test(NAME test_ff_pm_coupling_serial_implicit_stokes_first_stokes
              TARGET test_ff_reversed test_pm_reversed
              LABELS freeflow stokes precice darcy
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py 
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/stokes-iterative-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/stokes-iterative-00007.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/params.input 
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/precice-config.xml
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_coupling_serial_implicit_stokes_first_stokes PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_coupling_serial_implicit_stokes_first_darcy
              TARGET test_ff_reversed test_pm_reversed
              LABELS freeflow stokes precice darcy
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py 
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/darcy-iterative-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00007.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/params.input 
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/precice-config.xml
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_coupling_serial_implicit_stokes_first_darcy PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})


dumux_add_test(NAME test_ff_pm_coupling_serial_implicit_darcy_first_stokes
              TARGET test_ff_reversed test_pm_reversed
              LABELS freeflow stokes precice darcy
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py 
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/stokes-iterative-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/stokes-iterative-00007.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/params.input 
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/precice-config.xml
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_coupling_serial_implicit_darcy_first_stokes PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_coupling_serial_implicit_darcy_first_darcy
              TARGET test_ff_reversed test_pm_reversed
              LABELS freeflow stokes precice darcy
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py 
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/darcy-iterative-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00007.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/params.input 
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/precice-config.xml
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_coupling_serial_implicit_darcy_first_darcy PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_coupling_parallel_implicit_stokes
              TARGET test_ff_reversed test_pm_reversed
              LABELS freeflow stokes precice darcy
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py 
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/stokes-iterative-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/stokes-iterative-00008.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/params.input 
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/precice-config.xml
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_coupling_parallel_implicit_stokes PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_coupling_parallel_implicit_darcy
              TARGET test_ff_reversed test_pm_reversed
              LABELS freeflow stokes precice darcy
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py 
              CMD_ARGS --files ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/darcy-iterative-final.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/darcy-iterative-00008.vtu
              --dumux-param-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/params.input 
              --precice-config-file ${CMAKE_SOURCE_DIR}/test/coupling-ff-pm/reference-solutions/iterative/parallel-implicit/precice-config.xml
              --relative 1e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)

# Adding DuMuX directory to Python search path such that we can use DuMuX'x fuzzy (inexact) testing of vtu files
set_tests_properties(test_ff_pm_coupling_parallel_implicit_darcy PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/../dumux/bin/testing:$ENV{PYTHONPATH})


#CMD_ARGS       --script fuzzy
#                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel-reference.vtu
#                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_outflow-00002.vtu
#                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_channel params.input -Problem.OutletCondition Outflow
                             #-Problem.Name test_ff_stokes_channel_outflow"
