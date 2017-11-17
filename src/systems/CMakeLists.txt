set(systems_source_files
    system.cpp

    box.hpp
    system.hpp
)

add_library(systems STATIC ${systems_source_files})

target_link_libraries(systems particles enhance)