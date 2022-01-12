set(modules
"DumuxPreciceMacros.cmake"
"DumuxPreciceTestMacros.cmake"
)

include(GNUInstallDirs)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
