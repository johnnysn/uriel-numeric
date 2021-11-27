#include "pch.h"
#include <math.h> 
#include "method/ode/RushLarsenAdaptiveMethod.h"

#define ASSESS_LOCAL_ERROR \
	double f = (rhs[l] + rhs_new[l]) * 0.5; \
	double y_2nd_order = Y_old_[l] + dt * f; \
	double local_error = (abs(y_2nd_order) < EPSILON) ? abs(Y_new_[l] - EPSILON) : abs((y_2nd_order - Y_new_[l]) / y_2nd_order); \
	if (local_error > maxLocalError) maxLocalError = local_error


double RushLarsenAdaptiveMethod::step(
	double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double rel_tol, double dt_max
)
{
	CellModel* cellModel = (CellModel*)model;

	// a_{n}, b_{n} and f_{n} are (supposedly) already calculated from the previous step
	double* as_old = &(rhs[cellModel->HHStart]);
	double* bs_old = &(rhs[cellModel->nStates]);
	double* y_old_ = &(Y_old_[cellModel->HHStart]);
	double* y_new_ = &(Y_new_[cellModel->HHStart]);

	// Calculate the step y_{n+1}
	for (int l = 0; l < cellModel->nStates_HH; l++) {
		double a = as_old[l]; double b = bs_old[l];
		y_new_[l] = (abs(a) < EPSILON) ? (y_old_[l] + dt * (y_old_[l] * a + b)) : (exp(a * dt)*(y_old_[l] + b / a) - b / a);
	}
	for (int l = cellModel->MKStart; l < cellModel->MKEnd; l++)
		Y_new_[l] = Y_old_[l] + dt * rhs[l];
	for (int l = cellModel->NLStart; l < cellModel->NLEnd; l++)
		Y_new_[l] = Y_old_[l] + dt * rhs[l];

	// Calculate a_{n+1} and b_{n+1} and append them to the second half of the rhs array
	double* as = &(rhs[cellModel->HHStart + cellModel->nStates + cellModel->nStates_HH]);
	double* bs = &(rhs[cellModel->nStates + cellModel->nStates + cellModel->nStates_HH]);
	cellModel->calc_hh_coeff(as, bs, pars, algs, Y_new_, t + dt);
	// Calculate f_{n+1} and append it to the second half of the rhs array
	double* rhs_new = &(rhs[cellModel->nStates + cellModel->nStates_HH]);
	cellModel->calc_rhs_nl(rhs_new, pars, algs, Y_new_, t + dt);
	cellModel->calc_rhs_mk(rhs_new, pars, algs, Y_new_, t + dt);

	// Calculate max local error HH
	double maxLocalError = 0;
	for (int l = 0; l < cellModel->nStates_HH; l++) {
		double a = (as[l] + as_old[l]) * 0.5;
		double b = (bs[l] + bs_old[l]) * 0.5;
		double y_2nd_order = (abs(a) < EPSILON) ? y_old_[l] + dt * (y_old_[l] * a + b) : (exp(a * dt)*(y_old_[l] + b / a) - b / a);
		double local_error = (abs(y_2nd_order) < EPSILON) ? abs(y_new_[l] - EPSILON) : abs((y_2nd_order - y_new_[l]) / y_2nd_order);
		if (local_error > maxLocalError) maxLocalError = local_error;
	}
	for (int l = cellModel->NLStart; l < cellModel->NLEnd; l++) {
		ASSESS_LOCAL_ERROR;
	}
	for (int l = cellModel->MKStart; l < cellModel->MKEnd; l++) {
		ASSESS_LOCAL_ERROR;
	}
	
	return (maxLocalError > rel_tol) ? 
		-1 : (
				(maxLocalError < EPSILON) ? dt_max : std::fmin(dt * sqrt(0.5 * rel_tol / maxLocalError), dt_max)
			);
}

void RushLarsenAdaptiveMethod::prepare(Model* model, double* pars, double* algs, double* rhs, double* Y_ini_, double t_ini) {
	CellModel* cellModel = (CellModel*)model;

	cellModel->calc_rhs_nl(rhs, pars, algs, Y_ini_, t_ini);
	cellModel->calc_rhs_mk(rhs, pars, algs, Y_ini_, t_ini);

	double* as = &(rhs[cellModel->HHStart]);
	double* bs = &(rhs[cellModel->nStates]);
	cellModel->calc_hh_coeff(as, bs, pars, algs, Y_ini_, t_ini);
}

void RushLarsenAdaptiveMethod::updateRHS(Model* model, double* rhs) {
	CellModel* cellModel = (CellModel*)model;
	int separation_index = cellModel->nStates + cellModel->nStates_HH;
	for (int l = 0; l < separation_index; l++) rhs[l] = rhs[l + separation_index];
}