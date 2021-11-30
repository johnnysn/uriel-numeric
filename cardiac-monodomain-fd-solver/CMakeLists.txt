cmake_minimum_required (VERSION 2.6)
project (cardiac-monodomain-fd-solver)

include_directories("../uriel-numeric/include" "../uriel-tools/include")
set(SOURCES cardiac-monodomain-fd-solver.cpp )

add_definitions( -Wall -O3 -std=c++11 )
add_executable( cardiac-monodomain-fd-solver ${SOURCES} )
target_link_libraries(cardiac-monodomain-fd-solver ${CMAKE_SOURCE_DIR}/../uriel-numeric/release/liburiel-numeric.a ${CMAKE_SOURCE_DIR}/../uriel-tools/release/liburiel-tools.a)
target_compile_definitions(cardiac-monodomain-fd-solver PUBLIC LINUXOS)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
