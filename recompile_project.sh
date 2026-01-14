#!/bin/bash

function create_build_directory () {
	BASE_PATH=$1
	if [ -d "${BASE_PATH}/build" ]; then
	        echo "Found ${BASE_PATH}/build directory!"
	else
    		echo "Creating ${BASE_PATH}/build directory!"
        	mkdir -p "${BASE_PATH}/build"
	fi
}

function build_and_compile () {
	BASE_PATH=$1
	echo "$(tput setaf 7) $(tput setab 1) Building and compiling ${BASE_PATH} ... $(tput sgr 0)"
	cd ${BASE_PATH}/build
	cmake ..; make; cd ../..
}

# Build the 'uriel-numeric' and 'uriel-tools' static libraries with CMake
create_build_directory "uriel-numeric"
create_build_directory "uriel-tools"

build_and_compile "uriel-numeric"
build_and_compile "uriel-tools"

# Build the 'cardiac-cell-solver' with CMake
create_build_directory "cardiac-cell-solver"
build_and_compile "cardiac-cell-solver"

echo "$(tput setaf 7) $(tput setab 2) Sucessfull compilation! Binary for 'cardiac-cell-solver' generated at 'cardiac-cell-solver/bin'. $(tput sgr 0)"
