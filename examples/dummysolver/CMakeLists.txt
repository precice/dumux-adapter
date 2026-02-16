add_executable(dummy_participantOne participantOne.cc)
add_executable(dummy_participantTwo participantTwo.cc)

dune_symlink_to_source_files(FILES precice-dummy-solver-config.xml params_one.input params_two.input runFirst.sh runSecond.sh)
