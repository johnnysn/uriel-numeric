#pragma once

#include "model/ModelVariableDescription.h";

class Model {

private:

	// Array de descrições das variáveis
	ModelVariableDescription* descr;

public:

	Model(int nStates, int nAlgs, int nParams);

	// Number of states of the model
	int nStates;
	// Number of parameters of the model
	int nParams;
	// Number of algebraic variables of the model
	int nAlgs;

	virtual void calc_rhs(double* rhs, double* pars, double* algs, double* Y_old_, double t) = 0;

	void setStateDescr(int index, std::string name, std::string unit);
	
	std::string getStateName(int index);
	std::string getStateUnit(int index);

};