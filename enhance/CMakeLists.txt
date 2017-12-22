message(${CMAKE_CURRENT_SOURCE_DIR})

set(enhance_source_files
    logger.cpp
    random.cpp

    logger.hpp
    random.hpp
    
    singleton.hpp
    math_utility.hpp
    parallel.hpp
    observer_ptr.hpp
    triangular_matrix.hpp
)

add_library(enhance STATIC ${enhance_source_files})