#include "pch.h"
#include <math.h> 
#include "method/ode/RushLarsenMethod.h"

void RushLarsenMethod::step(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt) 
{
	CellModel* cellModel = (CellModel*)model;
	partitionedStep(Y_new_, model, pars, algs, rhs, Y_old_, t, dt);
	cellModel->calc_rhs_mk(rhs, pars, algs, Y_old_, t);
	for (int l = cellModel->MKStart; l < cellModel->MKEnd; l++)
		Y_new_[l] = Y_old_[l] + dt * rhs[l];
	cellModel->calc_rhs_nl(rhs, pars, algs, Y_old_, t);
	for (int l = cellModel->NLStart; l < cellModel->NLEnd; l++)
		Y_new_[l] = Y_old_[l] + dt * rhs[l];
}

void RushLarsenMethod::partitionedStep(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt) 
{
	CellModel* cellModel = (CellModel*)model;
	double* as = &(rhs[cellModel->HHStart]);
	double* bs = &(rhs[cellModel->nStates]);
	cellModel->calc_hh_coeff(as, bs, pars, algs, Y_old_, t);
	step(&(Y_new_[cellModel->HHStart]), cellModel->nStates_HH, as, bs, &(Y_old_[cellModel->HHStart]), dt);
}

void RushLarsenMethod::step(double* Y_new_, int n, double* as, double* bs, double* Y_old_, double dt) 
{
	for (int i = 0; i < n; i++) {
		if (as[i] == 0) { // TODO change to epsilon comparison
			Y_new_[i] = Y_old_[i] + dt * (Y_old_[i] * as[i] + bs[i]);
		} else {
			double aux = bs[i] / as[i];
			Y_new_[i] = exp(as[i] * dt)*(Y_old_[i] + aux) - aux;
		}
	}
}