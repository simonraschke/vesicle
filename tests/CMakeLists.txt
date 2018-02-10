find_package(Boost QUIET REQUIRED COMPONENTS system unit_test_framework)

add_executable(systems_test ${vesicle_SOURCE_DIR}/tests/systems_test.cpp)
target_link_libraries(systems_test systems IO enhance "-lboost_unit_test_framework")
add_test(NAME systems COMMAND systems_test)
target_include_directories(systems_test PUBLIC ${vesicle_SOURCE_DIR})
target_include_directories(systems_test PUBLIC ${vesicle_SOURCE_DIR}/src)

add_executable(particles_test ${vesicle_SOURCE_DIR}/tests/particles_test.cpp)
target_link_libraries(particles_test enhance "-lboost_unit_test_framework")
add_test(NAME particles COMMAND particles_test)
target_include_directories(particles_test PUBLIC ${vesicle_SOURCE_DIR})
target_include_directories(particles_test PUBLIC ${vesicle_SOURCE_DIR}/src)

add_executable(trajectory_reader_test ${vesicle_SOURCE_DIR}/tests/trajectory_reader_test.cpp)
target_link_libraries(trajectory_reader_test IO "-lboost_unit_test_framework")
add_test(NAME trajectory_reader COMMAND trajectory_reader_test)
target_include_directories(trajectory_reader_test PUBLIC ${vesicle_SOURCE_DIR})
target_include_directories(trajectory_reader_test PUBLIC ${vesicle_SOURCE_DIR}/src)
