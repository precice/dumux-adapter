set(modules
"DumuxPreciceMacros.cmake"
"DumuxPreciceTestMacros.cmake"
"AddPreciceConfigurationFileLinks.cmake"
)

include(GNUInstallDirs)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
