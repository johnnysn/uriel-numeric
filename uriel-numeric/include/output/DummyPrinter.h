#pragma once

#include "SolutionPrinter.h"

class DummyPrinter : public SolutionPrinter
{
public:
	void printNode(long nodeIndex, double time, int statec, double* statev);
	void printNodes(double time, int nodec, int nodeJump, int statec, double* statev);
};