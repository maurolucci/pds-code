
add_executable(converter converter.cpp)
target_link_libraries(converter PUBLIC mpgraphs pdslib)
enable_warnings(converter)