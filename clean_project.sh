#!/bin/bash

function remove_build_directory () {
	BASE_PATH=$1
	if [ -d "${BASE_PATH}/build" ]; then
	        echo "Removing ${BASE_PATH}/build directory!"
		rm -r ${BASE_PATH}/build
	else
    		echo "Cannot find ${BASE_PATH}/build directory!"
	fi
}

remove_build_directory "uriel-numeric"
remove_build_directory "uriel-tools"
remove_build_directory "cardiac-cell-solver"



