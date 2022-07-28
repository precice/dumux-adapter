install(FILES
	couplingadapter.hh
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux-precice)

add_library(dumux-precice STATIC couplingadapter.cc dumuxpreciceindexmapper.cc)
target_link_libraries(dumux-precice PRIVATE precice::precice)
