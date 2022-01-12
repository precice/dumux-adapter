add_executable(test_ff EXCLUDE_FROM_ALL main_ff.cc)
add_executable(test_pm EXCLUDE_FROM_ALL main_pm.cc)

target_compile_definitions(test_ff PUBLIC "ENABLEMONOLITHIC=0")
target_compile_definitions(test_pm PUBLIC "ENABLEMONOLITHIC=0")

add_input_file_links()
