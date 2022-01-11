install(FILES
	couplingadapter.hh
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux-precice)

add_library(dumux-precice SHARED couplingadapter.cc dumuxpreciceindexwrapper.cc)
target_link_libraries(dumux-precice PRIVATE precice::precice)

add_subdirectory(dumux-addon)
