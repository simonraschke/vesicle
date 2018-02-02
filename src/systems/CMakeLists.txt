message(${CMAKE_CURRENT_SOURCE_DIR})

find_package(Eigen3 REQUIRED NO_MODULE)

file(GLOB systems_source_files
    *.cpp
    *.hpp
)

add_library(systems STATIC ${systems_source_files})

target_link_libraries(systems particles enhance simulations thermostats Eigen3::Eigen)
target_include_directories(systems PUBLIC ${PROJECT_SOURCE_DIR}/src)