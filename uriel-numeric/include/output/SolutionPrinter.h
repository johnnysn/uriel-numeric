#pragma once

class SolutionPrinter
{
public:
	virtual void printNode(long nodeIndex, double time, int statec, double* statev) = 0;
	virtual void printNodes( double time, int nodec, int nodeJump, int statec, double* statev) = 0;
};