message(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(monte_carlo_utility)

file(GLOB simulations_source_files
    *.cpp
    *.hpp
)

add_library(simulations STATIC ${simulations_source_files})

target_link_libraries(simulations enhance mcutil)
target_include_directories(simulations PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(simulations PUBLIC ${PROJECT_SOURCE_DIR}/src)