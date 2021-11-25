#pragma once
#include "model/CellModel.h"

class RushLarsenMethod {

	/**
	* The Model must be an instance of CellModel, otherwise an error will be thrown
	*/
	void step(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt);
	void partitionedStep(double* Y_new_, Model* model, double* pars, double* algs, double* rhs, double* Y_old_, double t, double dt);

	void step(double* Y_new_, int n, double* as, double* bs, double* Y_old_, double dt);

private:
	CellModel* getCellModel(Model* model);

};