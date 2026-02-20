add_executable(dumuxprecice_filterInteriorTest filterInteriorTest.cc)
target_link_libraries(dumuxprecice_filterInteriorTest dunecommon dunegrid)

dumux_add_test(
  NAME dumuxprecice_filterInteriorTest
  TARGET dumuxprecice_filterInteriorTest
  LABELS unit grid
  MPI 2
)
