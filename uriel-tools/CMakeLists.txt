cmake_minimum_required (VERSION 2.6)
project (uriel-tools)

include_directories("include")
file(GLOB_RECURSE
	SOURCES 
	"source/pch.cpp"
	"source/options/*.cpp"
)
add_definitions( -Wall -O3 -std=c++11 )
ADD_LIBRARY( uriel-tools STATIC ${SOURCES} )
set_target_properties(uriel-tools
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "release"
    LIBRARY_OUTPUT_DIRECTORY "release"
    RUNTIME_OUTPUT_DIRECTORY "release"
)
target_compile_definitions(uriel-tools PUBLIC LINUXOS)
