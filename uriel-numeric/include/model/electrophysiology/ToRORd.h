#pragma once
#include "model/CellModel.h"

class ToRORd : public CellModel {
public:
	ToRORd();

	void calc_algs_nl(double* algs, double* pars, double* Y_old_, double time);
	void calc_algs_hh(double* algs, double* pars, double* Y_old_, double time);

	void calc_rhs_nl(double* rhs, double* pars, double* algs, double* Y_old_, double t);
	void calc_rhs_hh(double* rhs, double* pars, double* algs, double* Y_old_, double t);
	void calc_rhs_mk(double* rhs, double* pars, double* algs, double* Y_old_, double t);

	void calc_hh_coeff(double* a, double* b, double* pars, double* algs, double* Y_old_, double t);
	void prep_mk_transitions(double* algs, double* pars, double* Y_old_, double t);
	void calc_mk_transitions(double** T, int mk_index, double* pars, double* algs, double* Y_old_, double t);

	void set_default_parameters(double* pars);
	void set_default_initial_state(double* Y_old_);

	double calc_stimulus(double* pars, double t);
};