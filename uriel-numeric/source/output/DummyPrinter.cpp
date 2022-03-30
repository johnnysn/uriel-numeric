#include "pch.h"
#include "output/DummyPrinter.h"
#include <iostream>

void DummyPrinter::printNode(long nodeIndex, double time, int statec, double* statev)
{
	std::cout << "t = " << time << " node " << nodeIndex << ": y = [ ";

	for (int i = 0; i < statec; i++) {
		std::cout << statev[i] << " ";
	}

	std::cout << "]" << std::endl;
}

void DummyPrinter::printNodes(double time, int nodec, int nodeJump, int statec, double* statev)
{
	std::cout << "t = " << time << ",";
	
	for (int l = 0; l < nodec; l++) {
		int k = l * nodeJump;
		std::cout << " y = [ ";
		for (int i = k; i < k + statec; i++) {
			std::cout << statev[i] << " ";
		}
		std::cout << "]";
	}

	std::cout << std::endl;
}