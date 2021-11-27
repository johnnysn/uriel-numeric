#include "pch.h"
#include "method/ode/ForwardEulerMethod.h"

void ForwardEulerMethod::step(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt)
{
	model->calc_rhs(rhs, pars, algs, Y_old_, t);
	for (int l = 0; l < model->nStates; l++) {
		Y_new_[l] = Y_old_[l] + dt * rhs[l];
	}
}

void ForwardEulerMethod::partitionedStep(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt)
{
	step(Y_new_, model, pars, algs, rhs, Y_old_, t, dt);
}