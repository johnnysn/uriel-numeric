#include "pch.h"
#include <math.h> 
#include "method/ode/ForwardEulerAdaptiveMethod.h"

#define ASSESS_LOCAL_ERROR \
	double f = (rhs[l] + rhs_new[l]) * 0.5; \
	double y_2nd_order = Y_old_[l] + dt * f; \
	double local_error = (abs(y_2nd_order) < EPSILON) ? abs(Y_new_[l] - EPSILON) : abs((y_2nd_order - Y_new_[l]) / y_2nd_order); \
	if (local_error > maxLocalError) maxLocalError = local_error

double ForwardEulerAdaptiveMethod::step(
	double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double rel_tol, double dt_max
)
{
	// f_{n} is (supposedly) already calculated from the previous step
	// Calculate the step y_{n+1}
	for (int l = 0; l < model->nStates; l++)
		Y_new_[l] = Y_old_[l] + dt * rhs[l];

	// Calculate f_{n+1} and append it to the second half of the rhs array
	double* rhs_new = &(rhs[model->nStates]);
	model->calc_rhs(rhs_new, pars, algs, Y_new_, t + dt);

	// Calculate max local error HH
	double maxLocalError = 0;
	for (int l = 0; l < model->nStates; l++) {
		ASSESS_LOCAL_ERROR;
	}

	return (maxLocalError > rel_tol) ?
		-1 : (
			(maxLocalError < EPSILON) ? dt_max : std::fmin(dt * sqrt(0.5 * rel_tol / maxLocalError), dt_max)
			);
}

void ForwardEulerAdaptiveMethod::prepare(Model* model, double* pars, double* algs, double* rhs, double* Y_ini_, double t_ini) {
	model->calc_rhs(rhs, pars, algs, Y_ini_, t_ini);
}

void ForwardEulerAdaptiveMethod::updateRHS(Model* model, double* rhs) {
	for (int l = 0; l < model->nStates; l++) rhs[l] = rhs[l + model->nStates];
}