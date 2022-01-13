#pragma once

#include "model/Model.h"

class ODEMethod {

public:
	virtual void step(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double** strut = NULL) = 0;
	virtual void partitionedStep(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double** strut = NULL) = 0;

};