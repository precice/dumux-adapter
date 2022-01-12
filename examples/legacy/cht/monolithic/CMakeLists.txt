add_input_file_links()

add_executable(test_conjugate_heat_monolithic EXCLUDE_FROM_ALL main.cc)

target_compile_definitions(test_conjugate_heat_monolithic PUBLIC "ENABLEMONOLITHIC=1")
