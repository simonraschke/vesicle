message(${CMAKE_CURRENT_SOURCE_DIR})

set(systems_source_files
    box.cpp
    system.cpp
    controller.cpp
    
    box.hpp
    system.hpp
    controller.hpp
)

add_library(systems STATIC ${systems_source_files})

target_link_libraries(systems particles enhance simulations thermostats)
target_include_directories(systems PUBLIC ${PROJECT_SOURCE_DIR}/src)