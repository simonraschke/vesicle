add_executable(systems_test ${PROJECT_SOURCE_DIR}/tests/systems_test.cpp)
target_link_libraries(systems_test systems IO)
add_test(NAME systems COMMAND systems_test)
target_include_directories(systems_test PUBLIC ${PROJECT_SOURCE_DIR}/src)

add_executable(particles_test ${PROJECT_SOURCE_DIR}/tests/particles_test.cpp)
target_link_libraries(particles_test enhance)
add_test(NAME particles COMMAND particles_test)
target_include_directories(particles_test PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(particles_test PUBLIC ${PROJECT_SOURCE_DIR}/src)

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)c