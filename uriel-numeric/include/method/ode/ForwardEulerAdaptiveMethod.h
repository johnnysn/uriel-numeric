#pragma once
#include "model/Model.h"
#include "method/ode/ODEAdaptiveMethod.h"

class ForwardEulerAdaptiveMethod : public ODEAdaptiveMethod {

	double step(
		double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double rel_tol, double dt_max
	);

	void prepare(Model* model, double* pars, double* algs, double* rhs, double* Y_ini_, double t_ini);
	void updateRHS(Model* model, double* rhs);
private:

};