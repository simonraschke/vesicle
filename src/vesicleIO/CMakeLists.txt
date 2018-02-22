message(${CMAKE_CURRENT_SOURCE_DIR})

find_package(Boost COMPONENTS program_options system filesystem REQUIRED)

file(GLOB IO_source_files
    *.cpp
    *.hpp
)

add_library(IO STATIC ${IO_source_files})

# target_link_libraries(IO enhance Eigen::Eigen "-lboost_system -lboost_filesystem -lboost_program_options")
target_link_libraries(IO enhance "-lboost_system -lboost_filesystem -lboost_program_options")

target_include_directories(IO PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(IO PUBLIC ${PROJECT_SOURCE_DIR}/src)