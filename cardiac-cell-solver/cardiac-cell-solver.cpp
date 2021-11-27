// cardiac-cell-solver.cpp : Este arquivo contém a função 'main'. A execução do programa começa e termina ali.
//

#include <iostream>
#include "model/electrophysiology/Fox2002.h"
#include "model/electrophysiology/Tusscher2004.h"
#include "method/ode/RushLarsenAdaptiveMethod.h"
#include "method/ode/RushLarsenMethod.h"
#include "method/ode/ForwardEulerMethod.h"
#include "OptionParser.h"

#define METHOD_EULER 0
#define METHOD_RL 1
#define METHOD_AEULER 2
#define METHOD_ARL 3

using namespace std;

void solveFixed(ODEMethod* method, CellModel* model, double dt, double dt_save, double tf);
void solveADP(ODEAdaptiveMethod* method, CellModel* model, double dt, double dt_save, double tf, double dt_max, double rel_tol);

int main(int argc, char** argv)
{
	OptionParser::addOption("model", "Model: 0 -> ten Tusscher 2004, 1 -> Fox 2002");
	OptionParser::addOption("method", "Method: 0 -> Rush Larsen, 1 -> Rush Larsen ADP");
	OptionParser::addOption("dt", "Base time step.");
	OptionParser::addOption("dt_save", "Time step for saving.");
	OptionParser::addOption("tf", "Final time");
	OptionParser::addOption("dt_max", "Maximum time step for adaptive solvers.");
	OptionParser::addOption("rel_tol", "Relative tolerance for adaptive solvers.");

	OptionParser::parseOptions(argc, argv);

	int method_index = OptionParser::parseInt("method");
	int model_index = OptionParser::parseInt("model");

	double dt = OptionParser::foundOption("dt") ? OptionParser::parseDouble("dt") : 0.01;
	double dt_save = OptionParser::foundOption("dt_save") ? OptionParser::parseDouble("dt_save") : 1;
	double tf = OptionParser::foundOption("tf") ? OptionParser::parseDouble("tf") : 400;

	CellModel* model;
	if (model_index == 0) {
		model = new Tusscher2004();
	} else if (model_index == 1) {
		model = new Fox2002();
	} else {
		return -1;
	}

	if (method_index == METHOD_EULER) {
		solveFixed(new ForwardEulerMethod(), model, dt, dt_save, tf);
	} else if (method_index == METHOD_RL) {
		solveFixed(new RushLarsenMethod(), model, dt, dt_save, tf);
	} else if (method_index == METHOD_ARL) {
		double dt_max, rel_tol;
		if (OptionParser::foundOption("dt_max")) dt_max = OptionParser::parseDouble("dt_max");
		else dt_max = 0.1;

		if (OptionParser::foundOption("rel_tol")) rel_tol = OptionParser::parseDouble("rel_tol");
		else rel_tol = 0.02;

		solveADP(new RushLarsenAdaptiveMethod(), model, dt, dt_save, tf, dt_max, rel_tol);
	} else {
		return -1;
	}

    cout << "Simulation has ended.\n";
	return 0;
}

void solveFixed(ODEMethod* method, CellModel* model, double dt, double dt_save, double tf) {

	double* Y_old_ = new double[model->nStates];
	double* Y_new_ = new double[model->nStates];

	// rhs will store the righ-hand-side values of NL and MK ODEs, and the coefficients a and b of the HH equations 
	double* rhs = new double[model->nStates + model->nStates_HH];
	double* algs = new double[model->nAlgs];
	double* params = new double[model->nParams];

	model->set_default_initial_state(Y_old_);
	model->set_default_parameters(params);

	double t_save = 0; 
	for (double t = 0; t <= tf; t += dt) {
		method->step(Y_new_, model, params, algs, rhs, Y_old_, t, dt);

		for (int l = 0; l < model->nStates; l++) Y_old_[l] = Y_new_[l];

		t_save += dt;
		if (t_save >= dt_save) {
			cout << "t = " << t + dt << ", V = " << Y_new_[0] << endl;
			t_save = 0;
		}
	}
}

void solveADP(ODEAdaptiveMethod* method, CellModel* model, double dt, double dt_save, double tf, double dt_max, double rel_tol) {
	double* Y_old_ = new double[model->nStates];
	double* Y_new_ = new double[model->nStates];

	// rhs will store the righ-hand-side values of NL and MK ODEs, and the coefficients a and b of the HH equations 
	// For adaptive methods, the array must be doubled
	double* rhs = new double[2 * (model->nStates + model->nStates_HH)];
	double* algs = new double[model->nAlgs];
	double* params = new double[model->nParams];

	model->set_default_initial_state(Y_old_);
	model->set_default_parameters(params);
	method->prepare(model, params, algs, rhs, Y_old_, 0);

	double t_save = 0; double dt_new;
	for (double t = 0; t <= tf; t += dt, dt = dt_new) {
		dt_new = method->step(Y_new_, model, params, algs, rhs, Y_old_, t, dt, 0.02, 0.05);
		while (dt_new < 0) {
			dt = dt * 0.5;
			dt_new = method->step(Y_new_, model, params, algs, rhs, Y_old_, t, dt, 0.02, 0.05);
		}

		for (int l = 0; l < model->nStates; l++) Y_old_[l] = Y_new_[l];
		for (int l = 0; l < model->nStates + model->nStates_HH; l++) rhs[l] = rhs[l + model->nStates + model->nStates_HH];

		t_save += dt;
		if (t_save >= dt_save) {
			cout << "t = " << t + dt << ", V = " << Y_new_[0] << endl;
			t_save = 0;
		}
	}
}