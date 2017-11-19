message(${CMAKE_CURRENT_SOURCE_DIR})

set(simulations_source_files
    langevuin.cpp

    algorithm.hpp
    langevuin.hpp
)

add_library(simulations STATIC ${simulations_source_files})

target_link_libraries(simulations enhance)
target_include_directories(simulations PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(simulations PUBLIC ${PROJECT_SOURCE_DIR}/src)