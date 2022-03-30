#include "pch.h"
#include "output/SingleFilePrinter.h"
#include <iostream>

SingleFilePrinter::SingleFilePrinter(std::string fileName)
{
	fout.open(fileName);
}

SingleFilePrinter::~SingleFilePrinter()
{
	fout.close();
}

void SingleFilePrinter::printNode(long nodeIndex, double time, int statec, double* statev)
{
	fout << nodeIndex << " " << time << " ";

	for (int i = 0; i < statec; i++) {
		fout << statev[i] << " ";
	}

	fout << std::endl;
}

void SingleFilePrinter::printNodes(double time, int nodec, int nodeJump, int statec, double* statev)
{
	fout << time << " ";

	for (int l = 0; l < nodec; l++) {
		int k = l * nodeJump;
		for (int i = k; i < k + statec; i++) {
			fout << statev[i] << " ";
		}
	}
	fout << std::endl;
}
