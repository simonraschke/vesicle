message(${CMAKE_CURRENT_SOURCE_DIR})

set(enhance_source_files
    history.cpp
    logger.cpp
    random.cpp

    history.hpp
    logger.hpp
    random.hpp
    
    singleton.hpp
    math_utility.hpp
    parallel.hpp
    observer_ptr.hpp
)

add_library(enhance STATIC ${enhance_source_files})