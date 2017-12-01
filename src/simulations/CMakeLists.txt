message(${CMAKE_CURRENT_SOURCE_DIR})

set(simulations_source_files
    algorithm.cpp
    verlet.cpp
    shake_verlet.cpp
    langevin.cpp
    
    algorithm.hpp
    shake_verlet.hpp
    langevin.hpp
)

add_library(simulations STATIC ${simulations_source_files})

target_link_libraries(simulations enhance)
target_include_directories(simulations PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(simulations PUBLIC ${PROJECT_SOURCE_DIR}/src)