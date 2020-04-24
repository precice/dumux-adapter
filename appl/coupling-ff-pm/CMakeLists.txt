# add a target that builds all exercise solutions
add_custom_target(test_exercises)


add_subdirectory(monolithic)
add_subdirectory(fvca-monolithic-reversed)
add_subdirectory(iterative)
add_subdirectory(iterative-reversed)
add_subdirectory(iterative-reversed-mono-flux)
add_subdirectory(iterative-reversed-mono-flux-no-precice)
