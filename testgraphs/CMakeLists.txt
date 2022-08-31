
add_executable(converter converter.cpp)
target_link_libraries(converter PUBLIC setgraph pdslib)
enable_warnings(converter)