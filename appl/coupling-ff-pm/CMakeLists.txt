# add a target that builds all exercise solutions
add_custom_target(test_exercises)


add_subdirectory(monolithic-to-rename)
add_subdirectory(fvca-monolithic)
add_subdirectory(iterative)
add_subdirectory(fvca-iterative)
#add_subdirectory(iterative-reversed-mono-flux)
#add_subdirectory(iterative-reversed-mono-flux-no-precice)


#add_subdirectory(interface-coupling-test)

# monolithic
add_subdirectory(monolithic)
# partitioned
add_subdirectory(partitioned)