#include "pch.h"
#include "model/CellModel.h"

CellModel::CellModel(int nStates_NL, int nStates_HH, int nStates_MK, int nMarkovModels, int nAlgs, int nParams):
	Model(nStates_NL + nStates_HH + nStates_MK, nAlgs, nParams),
	nStates_NL(nStates_NL), nStates_HH(nStates_HH), nStates_MK(nStates_MK),
	nMarkovModels(nMarkovModels),
	NLStart(0), NLEnd(nStates_NL),
	HHStart(nStates_NL), HHEnd(nStates_NL + nStates_HH),
	MKStart(nStates_NL + nStates_HH), MKEnd(nStates_NL + nStates_HH + nStates_MK)
{
	nStates_MKM = (nStates_MK > 0) ? new int[nMarkovModels] : NULL;
}

void CellModel::calc_rhs(double* rhs, double* pars, double* algs, double* Y_old_, double t) 
{
	calc_rhs_nl(rhs, pars, algs, Y_old_, t);
	calc_rhs_hh(rhs, pars, algs, Y_old_, t);
	calc_rhs_mk(rhs, pars, algs, Y_old_, t);
}

