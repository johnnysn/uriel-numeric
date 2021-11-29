cmake_minimum_required (VERSION 2.6)
project (cardiac-cell-solver)

include_directories("../uriel-numeric/include" "../uriel-tools/include")
set(SOURCES cardiac-cell-solver.cpp )

add_definitions( -Wall -O3 -std=c++11 )
add_executable( cardiac-cell-solver ${SOURCES} )
target_link_libraries(cardiac-cell-solver ${CMAKE_SOURCE_DIR}/../uriel-numeric/release/liburiel-numeric.a ${CMAKE_SOURCE_DIR}/../uriel-tools/release/liburiel-tools.a)
target_compile_definitions(cardiac-cell-solver PUBLIC LINUXOS)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
