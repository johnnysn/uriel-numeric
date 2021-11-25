#pragma once

#include "model/Model.h"

class ODEAdaptiveMethod 
{

public:
	virtual void prepare(Model* model, double* pars, double* algs, double* rhs, double* Y_ini_, double t_ini) = 0;
	virtual double step(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double rel_tol, double dt_max) = 0;

};
