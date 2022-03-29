#pragma once

#include "SolutionPrinter.h"

class DummyPrinter : public SolutionPrinter
{
public:
	void printNode(long nodeIndex, double time, int statec, double* statev);
};