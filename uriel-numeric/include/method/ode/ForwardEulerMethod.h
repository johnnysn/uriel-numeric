#pragma once
#include "model/Model.h"
#include "method/ode/ODEMethod.h"

class ForwardEulerMethod : public ODEMethod {

	void step(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double** strut = NULL);
	void partitionedStep(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double** strut = NULL);

};