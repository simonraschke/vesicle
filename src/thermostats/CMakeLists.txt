message(${CMAKE_CURRENT_SOURCE_DIR})

set(thermostats_source_files
    thermostat.cpp
    andersen.cpp
    
    thermostat.hpp
    andersen.hpp
)

add_library(thermostats STATIC ${thermostats_source_files})

# target_link_libraries(thermostats )
target_include_directories(thermostats PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(thermostats PUBLIC ${PROJECT_SOURCE_DIR}/src)