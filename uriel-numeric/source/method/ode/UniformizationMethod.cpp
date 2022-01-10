#include "pch.h"
#include <math.h> 
#include "method/ode/UniformizationMethod.h"

UniformizationMethod::UniformizationMethod() {
	rushLarsenMethod = new RushLarsenMethod();
}

UniformizationMethod::~UniformizationMethod() {
	delete rushLarsenMethod;
}

void UniformizationMethod::step(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt)
{
	CellModel* cellModel = (CellModel*)model;
	rushLarsenMethod->partitionedStep(Y_new_, model, pars, algs, rhs, Y_old_, t, dt);
	partitionedStep(Y_new_, model, pars, algs, rhs, Y_old_, t, dt);
	cellModel->calc_rhs_nl(rhs, pars, algs, Y_old_, t);
	for (int l = cellModel->NLStart; l < cellModel->NLEnd; l++)
		Y_new_[l] = Y_old_[l] + dt * rhs[l];
}

void UniformizationMethod::partitionedStep(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt)
{
	CellModel* cellModel = (CellModel*)model;
	int startIndex = cellModel->MKStart;

	cellModel->prep_mk_transitions(algs, pars, Y_old_, t);
	for (int k = 0; k < cellModel->nMarkovModels; k++) {
		double* y_old_ = &(Y_old_[startIndex]);
		double* y_new_ = &(Y_new_[startIndex]);
		int transitionMatrixIndex = cellModel->nStates + cellModel->nStates_HH;
		double* T = &(rhs[transitionMatrixIndex]);
		cellModel->calc_mk_transitions(T, k, pars, algs, Y_old_, t);

		step(y_new_, cellModel->nStates_MKM[k], T, y_old_, dt);

		startIndex += cellModel->nStates_MKM[k];
	}
}

void UniformizationMethod::step(double* Y_new_, int n, double* T, double* Y_old_, double dt) {

}