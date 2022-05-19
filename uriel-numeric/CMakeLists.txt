cmake_minimum_required (VERSION 2.6)
project (uriel-numeric)

include_directories("include")
file(GLOB_RECURSE
	SOURCES 
	"source/pch.cpp"
	"source/output/*.cpp"
	"source/method/ode/*.cpp"
	"source/model/*.cpp"
	"source/model/electrophysiology/*.cpp"
)
add_definitions( -Wall -O3 -std=c++11 )
ADD_LIBRARY( uriel-numeric STATIC ${SOURCES} )
set_target_properties(uriel-numeric
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "release"
    LIBRARY_OUTPUT_DIRECTORY "release"
    RUNTIME_OUTPUT_DIRECTORY "release"
)
target_compile_definitions(uriel-numeric PUBLIC LINUXOS)
