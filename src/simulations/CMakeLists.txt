message(${CMAKE_CURRENT_SOURCE_DIR})

set(simulations_source_files
    algorithm.cpp
    verlet.cpp
    shake_verlet.cpp
    langevin.cpp
    monte_carlo.cpp
    
    algorithm.hpp
    verlet.hpp
    shake_verlet.hpp
    langevin.hpp
    monte_carlo.hpp
)

add_library(simulations STATIC ${simulations_source_files})

target_link_libraries(simulations enhance)
target_include_directories(simulations PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(simulations PUBLIC ${PROJECT_SOURCE_DIR}/src)