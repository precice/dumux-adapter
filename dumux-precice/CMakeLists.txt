install(FILES
	couplingadapter.hh
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux-precice)

dune_library_add_sources(dumux-precice SOURCES couplingadapter.cc dumuxpreciceindexmapper.cc)
target_link_libraries(dumux-precice PRIVATE precice::precice)
