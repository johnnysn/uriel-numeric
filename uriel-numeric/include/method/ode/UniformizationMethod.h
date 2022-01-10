#pragma once
#include "model/CellModel.h"
#include "method/ode/ODEMethod.h"
#include "method/ode/RushLarsenMethod.h"

class UniformizationMethod : public ODEMethod {

public:

	UniformizationMethod();
	virtual ~UniformizationMethod();

	/**
	* The Model must be an instance of CellModel, otherwise an error will be thrown
	*/
	void step(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt);
	void partitionedStep(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt);

private:

	RushLarsenMethod* rushLarsenMethod;
	void step(double* Y_new_, int n, double* T, double* Y_old_, double dt);

};