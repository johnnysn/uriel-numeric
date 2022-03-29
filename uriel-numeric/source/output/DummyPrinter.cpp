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