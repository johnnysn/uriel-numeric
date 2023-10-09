// cardiac-cell-solver.cpp : Este arquivo contém a função 'main'. A execução do programa começa e termina ali.
//

#include <iostream>
#include "model/electrophysiology/Fox2002.h"
#include "model/electrophysiology/Tusscher2004.h"
#include "model/electrophysiology/Bondarenko2004.h"
#include "model/electrophysiology/Noble1962.h"
#include "model/electrophysiology/ToRORd_fkatp_2019.h"
#include "method/ode/RushLarsenAdaptiveMethod.h"
#include "method/ode/ForwardEulerAdaptiveMethod.h"
#include "method/ode/RushLarsenMethod.h"
#include "method/ode/UniformizationMethod.h"
#include "method/ode/ForwardEulerMethod.h"
#include "output/DummyPrinter.h"
#include "output/SingleFilePrinter.h"
#include "options/OptionParser.h"
#include "options/IndexedValue.h"

#define METHOD_EULER 0
#define METHOD_RL 1
#define METHOD_AEULER 2
#define METHOD_ARL 3
#define METHOD_UNI 4

using namespace std;

void solveFixed(ODEMethod* method, CellModel* model, double dt, double dt_save, double tf, SolutionPrinter* printer, vector<IndexedValue> parameters);
void solveADP(ODEAdaptiveMethod* method, CellModel* model, double dt, double dt_save, SolutionPrinter* printer, double tf, double dt_max, double rel_tol, vector<IndexedValue> parameters);

int main(int argc, char** argv)
{
	OptionParser::addOption("model", "Model: 0 -> ten Tusscher 2004, 1 -> Fox 2002, 2 -> Bondarenko 2004");
	OptionParser::addOption("method", "Method: 0 -> Euler, 1 -> Rush Larsen, 2 -> Euler ADP, 3 -> Rush Larsen ADP, 4 -> UNI");
	OptionParser::addOption("dt", "Base time step.");
	OptionParser::addOption("dt_save", "Time step for saving.");
	OptionParser::addOption("tf", "Final time");
	OptionParser::addOption("dt_max", "Maximum time step for adaptive solvers.");
	OptionParser::addOption("rel_tol", "Relative tolerance for adaptive solvers.");
	OptionParser::addOption("outputFile", "Filename for printing output");
	OptionParser::addOption("parameters", "Set extra parameters");

	OptionParser::parseOptions(argc, argv);

	int method_index = OptionParser::parseInt("method");
	int model_index = OptionParser::parseInt("model");

	double dt = OptionParser::foundOption("dt") ? OptionParser::parseDouble("dt") : 0.01;
	double dt_save = OptionParser::foundOption("dt_save") ? OptionParser::parseDouble("dt_save") : 1;
	double tf = OptionParser::foundOption("tf") ? OptionParser::parseDouble("tf") : 400;

	vector<IndexedValue> parameters;
	if (OptionParser::foundOption("parameters")) {
		string parameters_str = OptionParser::optionValue("parameters");
		parameters = OptionParser::parseIndexedValues("parameters");
		cout << parameters_str << endl;
		for (unsigned int i = 0; i < parameters.size(); i++) {
			cout << parameters[i].index << ": " << parameters[i].value << endl;
		}
	}

	CellModel* model;
	if (model_index == 0) {
		model = new Tusscher2004();
	} else if (model_index == 1) {
		model = new Fox2002();
	} else if (model_index == 2) {
		model = new Bondarenko2004();
	} else {
		cout << "Invalid model index." << endl;
		return EXIT_FAILURE;
	}

	SolutionPrinter* printer;
	if (OptionParser::foundOption("outputFile"))
		printer = new SingleFilePrinter(OptionParser::optionValue("outputFile"));
	else
		printer = new DummyPrinter();

	if (method_index == METHOD_EULER) {
		solveFixed(new ForwardEulerMethod(), model, dt, dt_save, tf, printer, parameters);
	} else if (method_index == METHOD_RL) {
		solveFixed(new RushLarsenMethod(), model, dt, dt_save, tf, printer, parameters);
	} else if (method_index == METHOD_UNI) {
		solveFixed(new UniformizationMethod(), model, dt, dt_save, tf, printer, parameters);
	} else if (method_index == METHOD_ARL || method_index == METHOD_AEULER) {
		double dt_max, rel_tol;
		if (OptionParser::foundOption("dt_max")) dt_max = OptionParser::parseDouble("dt_max");
		else dt_max = 0.1;

		if (OptionParser::foundOption("rel_tol")) rel_tol = OptionParser::parseDouble("rel_tol");
		else rel_tol = 0.02;

		ODEAdaptiveMethod* method;
		if (method_index == METHOD_ARL) method = new RushLarsenAdaptiveMethod();
		else method = new ForwardEulerAdaptiveMethod();

		solveADP(method, model, dt, dt_save, printer, tf, rel_tol, dt_max, parameters);
	} else {
		cout << "Invalid method index." << endl;
		return EXIT_FAILURE;
	}

	delete printer;
    cout << "Simulation has ended.\n";
	return 0;
}

void solveFixed(ODEMethod* method, CellModel* model, double dt, double dt_save, double tf, SolutionPrinter* printer, vector<IndexedValue> parameters)
{

	double* Y_old_ = new double[model->nStates];
	double* Y_new_ = new double[model->nStates];

	double** Tr = NULL;
	if (model->nStates_MKM_max > 0) {
		Tr = new double* [model->nStates_MKM_max];
		for (int i = 0; i < model->nStates_MKM_max; ++i)
			Tr[i] = new double[model->nStates_MKM_max];
	}

	// rhs will store the righ-hand-side values of NL and MK ODEs, and the coefficients a and b of the HH equations (for the RL method)
	double* rhs = new double[model->nStates + model->nStates_HH];
	double* algs = new double[model->nAlgs];
	double* params = new double[model->nParams];


	model->set_default_initial_state(Y_old_);
	model->set_default_parameters(params);
	for (unsigned int i = 0; i < parameters.size(); i++) {
		params[parameters[i].index] = parameters[i].value;
	}

	double t_save = 0; 
	for (double t = 0; t <= tf; t += dt) {
		method->step(Y_new_, model, params, algs, rhs, Y_old_, t, dt, Tr);

		for (int l = 0; l < model->nStates; l++) Y_old_[l] = Y_new_[l];

		t_save += dt;
		if (t_save >= dt_save) {
			printer->printNode(0, t + dt, 1, Y_new_);
			t_save = 0;
		}
	}
}

void solveADP(ODEAdaptiveMethod* method, CellModel* model, double dt, double dt_save, SolutionPrinter* printer, double tf, double rel_tol, double dt_max, vector<IndexedValue> parameters)
{
	double* Y_old_ = new double[model->nStates];
	double* Y_new_ = new double[model->nStates];

	// rhs will store the righ-hand-side values of NL and MK ODEs, and the coefficients a and b of the HH equations 
	// For adaptive methods, the array must be doubled
	double* rhs = new double[2 * (model->nStates + model->nStates_HH)];
	double* algs = new double[model->nAlgs];
	double* params = new double[model->nParams];

	model->set_default_initial_state(Y_old_);
	model->set_default_parameters(params);
	for (unsigned int i = 0; i < parameters.size(); i++) {
		params[parameters[i].index] = parameters[i].value;
	}

	method->prepare(model, params, algs, rhs, Y_old_, 0);

	double t_save = 0; double dt_new;
	for (double t = 0; t <= tf; t += dt, dt = dt_new) {
		dt_new = method->step(Y_new_, model, params, algs, rhs, Y_old_, t, dt, rel_tol, dt_max);
		while (dt_new < 0) {
			dt = dt * 0.5;
			dt_new = method->step(Y_new_, model, params, algs, rhs, Y_old_, t, dt, rel_tol, dt_max);
		}

		for (int l = 0; l < model->nStates; l++) Y_old_[l] = Y_new_[l];
		method->updateRHS(model, rhs);

		t_save += dt;
		if (t_save >= dt_save) {
			printer->printNode(0, t + dt, 1, Y_new_);
			t_save = 0;
		}
	}
}