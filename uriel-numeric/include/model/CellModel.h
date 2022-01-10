#pragma once

#include "model/Model.h"

/*
	Computational representation of a mathematical cell model. The variables must be orderly separated in three types:
	1) Fully nonlinear states
	2) Hodgkin-Huxley-type states
	3) Markov-type states
*/
class CellModel: public Model {

public:
	// Number of states described by fully nonlinear DEs
	const int nStates_NL;

	// Number of states described Hodgkin-Huxley-type ODEs
	const int nStates_HH;

	// Number of Markov-based states
	const int nStates_MK;

	// Number of Markov models
	const int nMarkovModels;

	// Number of states of each Markov model
	int* nStates_MKM;

	CellModel(int nStates_NL, int nStates_HH, int nStates_MK, int nMarkovModels, int nAlgs, int nParams);

	// Implemented inherited methods
	void calc_rhs(double* rhs, double* pars, double* algs, double* Y_old_, double t);

	// Abstract methods
	virtual void calc_rhs_nl(double* rhs, double* pars, double* algs, double* Y_old_, double t) = 0;
	virtual void calc_rhs_hh(double* rhs, double* pars, double* algs, double* Y_old_, double t) = 0;
	virtual void calc_rhs_mk(double* rhs, double* pars, double* algs, double* Y_old_, double t) = 0;

	virtual void calc_hh_coeff(double* a, double* b, double* pars, double* algs, double* Y_old_, double t) = 0;
	virtual void prep_mk_transitions(double* algs, double* pars, double* Y_old_, double t) = 0;
	virtual void calc_mk_transitions(double* T, int mk_index, double* pars, double* algs, double* Y_old_, double t) = 0;

	virtual void set_default_parameters(double* pars) = 0;
	virtual void set_default_initial_state(double* Y_old_) = 0;

	const int NLStart; const int NLEnd;
	const int HHStart; const int HHEnd;
	const int MKStart; const int MKEnd;

};