// cardiac-monodomain-fd-solver.cpp : Este arquivo contém a função 'main'. A execução do programa começa e termina ali.
//

#include <iostream>
#include "model/electrophysiology/Fox2002.h"
#include "model/electrophysiology/Tusscher2004.h"
#include "model/electrophysiology/Bondarenko2004.h"
#include "method/ode/RushLarsenMethod.h"
#include "method/ode/UniformizationMethod.h"
#include "options/OptionParser.h"
#include "output/DummyPrinter.h"
#include "output/SingleFilePrinter.h"
#include <chrono>
#include <omp.h>

long Nx, Ny;
CellModel* model;

long offset(long xi, long yi)
{
	if (yi == -1) yi = 1;
	else if (yi == Ny) yi = Ny - 2;
	if (xi == -1) xi = 1;
	else if (xi == Nx) xi = Nx - 2;
	return (Nx * yi + xi) * model->nStates;
}

void deleteStructures(double* Y_old_, double* Y_new_, double* ALGS, double* RHS, double* PARAMS, double** Tr)
{
	delete[] Y_old_;
	delete[] Y_new_;
	delete[] ALGS;
	delete[] RHS;
	delete[] PARAMS;
	if (Tr != NULL) {
		for (int i = 0; i < Nx * Ny * model->nStates_MKM_max; ++i) delete[] Tr[i];
		delete[] Tr;
	}
}


int main(int argc, char** argv)
{
	OptionParser::addOption("model", "Model: 0 -> ten Tusscher 2004, 1 -> Fox 2002, 2 -> Bondarenko 2004");
	OptionParser::addOption("method", "Method: 0 -> Rush Larsen Method, 1 -> Uniformization Method");
	OptionParser::addOption("dt", "Base time step.");
	OptionParser::addOption("dt_save", "Time step for saving.");
	OptionParser::addOption("tf", "Final time");
	OptionParser::addOption("dx", "Spatial discretization step in the x-direction (cm)");
	OptionParser::addOption("dy", "Spatial discretization step in the y-direction (cm)");
	OptionParser::addOption("Lx", "Horizontal length in cm");
	OptionParser::addOption("Ly", "Vertical length in cm");
	OptionParser::addOption("threads", "Number of threads");
	OptionParser::addOption("outputFile", "Filename for printing output");

	OptionParser::parseOptions(argc, argv);

	int model_index = OptionParser::parseInt("model");
	int method_index = OptionParser::parseInt("method");
	int num_threads = OptionParser::foundOption("threads") ? OptionParser::parseInt("threads") : 1;
	double dt = OptionParser::foundOption("dt") ? OptionParser::parseDouble("dt") : 0.02;
	double dt_save = OptionParser::foundOption("dt_save") ? OptionParser::parseDouble("dt_save") : 1;
	double tf = OptionParser::foundOption("tf") ? OptionParser::parseDouble("tf") : 400;
	double dx = OptionParser::foundOption("dx") ? OptionParser::parseDouble("dx") : 0.02;
	double dy = OptionParser::foundOption("dy") ? OptionParser::parseDouble("dy") : 0.02;
	double Lx = OptionParser::foundOption("Lx") ? OptionParser::parseDouble("Lx") : 2;
	double Ly = OptionParser::foundOption("Ly") ? OptionParser::parseDouble("Ly") : 2;

	omp_set_num_threads(num_threads);

	// Monodomain params
	double sigma_x = 1.0 / 163.0, sigma_y = 1.0 / 163.0, chi = 4, Cm = 1;

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

	ODEMethod* method;
	if (method_index == 0) {
		method = new RushLarsenMethod();
	} else if (method_index == 1) {
		method = new UniformizationMethod();
	} else {
		cout << "Invalid method index." << endl;
		return EXIT_FAILURE;
	}

	Nx = Lx / dx;
	Ny = Ly / dy;

	double* Y_old_ = new double [Nx * Ny * model->nStates]; 
	double* Y_new_ = new double [Nx * Ny * model->nStates];
	double* ALGS = new double [Nx * Ny * model->nAlgs]; double* PARAMS = new double [Nx * Ny * model->nParams];
	double* RHS = new double [Nx * Ny * (model->nStates + model->nStates_HH)];

	double** Tr = NULL;
	if (model->nStates_MKM_max > 0) {
		Tr = new double* [Nx * Ny * model->nStates_MKM_max];
		for (int i = 0; i < Nx * Ny * model->nStates_MKM_max; ++i)
			Tr[i] = new double[model->nStates_MKM_max];
	}

	double t = 0, t_save = 0;

	SolutionPrinter* printer;
	if (OptionParser::foundOption("outputFile"))
		printer = new SingleFilePrinter(OptionParser::optionValue("outputFile"));
	else
		printer = new DummyPrinter();

	auto start = std::chrono::high_resolution_clock::now();

	#pragma omp parallel 
	{

		#pragma omp for
		for (int k = 0; k < Nx * Ny; k++) {
			double* params = &(PARAMS[k * model->nParams]);
			model->set_default_parameters(params);
			params[0] = -1;// stim_state set to off
			model->set_default_initial_state(&(Y_old_[k * model->nStates]));
		}

		while (t < tf) {
			#pragma omp for
			for (long yi = 0; yi < Ny; yi++) {
				for (long xi = 0; xi < Nx; xi++) {
					long k = yi * Nx + xi;
					long m = k * model->nStates;
					double* y_new_ = &(Y_new_[m]); double* y_old_ = &(Y_old_[m]);
					double* algs = &(ALGS[k * model->nAlgs]); double* params = &(PARAMS[k * model->nParams]);
					double* rhs = &(RHS[k * (model->nStates + model->nStates_HH)]);
					double** tr = Tr != NULL ? &(Tr[k * model->nStates_MKM_max]) : NULL;

					// Stimulus
					double x = xi * dx, y = yi * dy;
					params[0] = (t > 1 && t < 3 && x > 0 && x < 0.1) ? 1 : -1;

					method->step(y_new_, model, params, algs, rhs, y_old_, t, dt, tr);

					//Laplacian term
					y_new_[0] += sigma_x / (Cm * chi) * dt / (dx * dx) *
						(Y_old_[offset(xi - 1, yi)] - 2 * y_old_[0] + Y_old_[offset(xi + 1, yi)]) +
						sigma_y / (Cm * chi) * dt / (dy * dy) *
						(Y_old_[offset(xi, yi - 1)] - 2 * y_old_[0] + Y_old_[offset(xi, yi + 1)]);
				}
			}

			#pragma omp single
			{
				double* aux = Y_old_;
				Y_old_ = Y_new_;
				Y_new_ = aux;

				t += dt;
				t_save += dt;
				if (t_save >= dt_save) {
					printer->printNodes(t, 2, (Nx*Ny - 1)*model->nStates, 1, Y_old_);
					t_save = 0;
				}
			}
		}
	}

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;

	delete printer;
	deleteStructures(Y_old_, Y_new_, ALGS, RHS, PARAMS, Tr);

    std::cout << "Simulation finished in " << elapsed.count() << " s." << std::endl;
}