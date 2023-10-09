#pragma once

#include "SolutionPrinter.h"
#include <iostream>
#include <fstream>

class SingleFilePrinter : public SolutionPrinter
{
public:
	SingleFilePrinter(std::string fileName);
	~SingleFilePrinter();

	void printNode(long nodeIndex, double time, int statec, double* statev);
	void printNodes(double time, int nodec, int nodeJump, int statec, double* statev);
private:
	std::ofstream fout;
};