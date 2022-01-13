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
	void step(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double** Tr);
	void partitionedStep(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double** Tr);

private:

	RushLarsenMethod* rushLarsenMethod;
	void step(double* Y_new_, int n, double* aux, double** Tr, double* pi, double* Y_old_, double dt);

	int calcNMax(double q, double dt);

	double tol = 1e-6;

};