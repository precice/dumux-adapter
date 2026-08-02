add_executable(dumuxprecice_filterInteriorTest EXCLUDE_FROM_ALL filterInteriorTest.cc)

dumux_add_test(
  NAME dumuxprecice_filterInteriorTest
  TARGET dumuxprecice_filterInteriorTest
  LABELS unit grid
  MPI 2
)
