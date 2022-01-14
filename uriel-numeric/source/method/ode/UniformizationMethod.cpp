#include "pch.h"
#include <math.h> 
#include "method/ode/UniformizationMethod.h"

UniformizationMethod::UniformizationMethod() {
	rushLarsenMethod = new RushLarsenMethod();
}

UniformizationMethod::~UniformizationMethod() {
	delete rushLarsenMethod;
}

void UniformizationMethod::step(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double** Tr)
{
	CellModel* cellModel = (CellModel*)model;
	cellModel->calc_rhs_nl(rhs, pars, algs, Y_old_, t);
	for (int l = cellModel->NLStart; l < cellModel->NLEnd; l++)
		Y_new_[l] = Y_old_[l] + dt * rhs[l];
	rushLarsenMethod->partitionedStep(Y_new_, model, pars, algs, rhs, Y_old_, t, dt);
	partitionedStep(Y_new_, model, pars, algs, rhs, Y_old_, t, dt, Tr);
}

void UniformizationMethod::partitionedStep(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt, double** Tr)
{
	CellModel* cellModel = (CellModel*)model;
	int startIndex = cellModel->MKStart;

	cellModel->prep_mk_transitions(algs, pars, Y_old_, t);

	double* aux = new double[cellModel->nStates_MKM_max];
	for (int k = 0; k < cellModel->nMarkovModels; k++) {
		double* y_old_ = &(Y_old_[startIndex]);
		double* y_new_ = &(Y_new_[startIndex]);
		double* pi = &(rhs[startIndex]);
		int m = cellModel->nStates_MKM[k];

		for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) Tr[i][j] = 0;
		cellModel->calc_mk_transitions(Tr, k, pars, algs, Y_old_, t);

		step(y_new_, m, aux, Tr, pi, y_old_, dt);

		startIndex += m;
	}
	delete[] aux;
}

void UniformizationMethod::step(double* Y_new_, int m, double* aux, double** Tr, double* pi, double* Y_old_, double dt)
{
	// Calculate q
	double q = 0;
	for (int j = 0; j < m; j++) {
		double sum = 0;
		Tr[j][j] = 0;
		for (int i = 0; i < m; i++) if (i != j) {
			sum += abs(Tr[i][j]);
			Tr[j][j] -= Tr[i][j];
		}
		if (sum > q) q = sum;
	}
	q += 1.0;

	// Calculate matrix Q and prepare aux array
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			Tr[i][j] = Tr[i][j] / q;
			if (i == j) Tr[i][j] += 1.0;
		}
		Y_new_[i] = aux[i] = Y_old_[i];
	}

	// Y_new_ is already Y_old_, so the sum starts with k=1
	double c = exp(-q * dt);
	double coeff = 1.0, coeff_sum = 1.0;
	for (int k = 1; k <= N_MAX; k++) {
		coeff *= q * dt / (double)k;

		// pi = Tr * aux
		for (int i = 0; i < m; i++) {
			pi[i] = 0;
			for (int j = 0; j < m; j++)
				pi[i] += Tr[i][j] * aux[j];
		}

		for (int i = 0; i < m; i++) {
			Y_new_[i] += coeff * pi[i];
			aux[i] = pi[i];
		}

		coeff_sum += coeff;
		if (tol > (1.0 - c * coeff_sum)) break; // Maximum N for tol is reached
	}

	// Multiply Z by exp(-q * dt)
	for (int i = 0; i < m; i++) {
		Y_new_[i] = c * Y_new_[i];
	}
}