#pragma once
#include "model/CellModel.h"
#include "method/ode/ODEAdaptiveMethod.h"

class RushLarsenAdaptiveMethod : public ODEAdaptiveMethod {

	/**
	* The Model must be an instance of CellModel, otherwise an error will be thrown
	*/
	double step(
		double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double rel_tol, double dt_max
	);

	void prepare(Model* model, double* pars, double* algs, double* rhs, double* Y_ini_, double t_ini);
private:

};