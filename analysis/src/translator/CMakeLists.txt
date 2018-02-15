message(${CMAKE_CURRENT_SOURCE_DIR})

# find_package(Eigen3 REQUIRED NO_MODULE)

file(GLOB translator_source_files
    *.cpp
    *.hpp
)

add_library(translator STATIC ${translator_source_files})
set_target_properties(translator PROPERTIES CXX_STANARD 14 )

target_link_libraries(translator particles)
target_include_directories(translator PUBLIC ${PROJECT_SOURCE_DIR}/src)