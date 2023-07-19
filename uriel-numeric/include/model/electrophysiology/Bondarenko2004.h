#pragma once
#include "model/CellModel.h"

class Bondarenko2004 : public CellModel {
public:
	Bondarenko2004();

	void calc_algs_nl(double* algs, double* pars, double* Y_old_, double time);
	void calc_algs_hh(double* algs, double* pars, double* Y_old_, double time);

	void calc_rhs_nl(double* rhs, double* pars, double* algs, double* Y_old_, double t);
	void calc_rhs_hh(double* rhs, double* pars, double* algs, double* Y_old_, double t);
	void calc_rhs_mk(double* rhs, double* pars, double* algs, double* Y_old_, double t, int mk_index = -1);

	void calc_hh_coeff(double* a, double* b, double* pars, double* algs, double* Y_old_, double t);
	void prep_mk_transitions(double* algs, double* pars, double* Y_old_, double t, int mk_index = -1);
	void calc_mk_transitions(double** T, int mk_index, double* pars, double* algs, double* Y_old_, double t);

	void set_default_parameters(double* pars);
	void set_default_initial_state(double* Y_old_);

	double calc_stimulus(double* pars, double t);

	bool has_single_rhs_formula(int i);
	double calc_single_rhs_formula(int i, double* pars, double* Y_old_, double t);
	bool is_mk_active(int i);
};