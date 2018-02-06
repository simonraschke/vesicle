message(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB enhance_source_files
    *.cpp
    *.hpp
)

add_library(enhance STATIC ${enhance_source_files})

target_compile_options(enhance PRIVATE -std=c++14 -O3 -g0 -Weffc++ -Wall -Wextra -Wpedantic -march=nocona -mtune=generic)
target_compile_definitions(enhance PUBLIC NDEBUG EIGEN_NO_DEBUG)