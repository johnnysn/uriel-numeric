// cardiac-monodomain-fd-solver.cpp : Este arquivo contém a função 'main'. A execução do programa começa e termina ali.
//

#include <iostream>
#include "model/electrophysiology/Fox2002.h"
#include "model/electrophysiology/Tusscher2004.h"
#include "method/ode/RushLarsenMethod.h"
#include "options/OptionParser.h"
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

int main(int argc, char** argv)
{
	OptionParser::addOption("model", "Model: 0 -> ten Tusscher 2004, 1 -> Fox 2002");
	OptionParser::addOption("dt", "Base time step.");
	OptionParser::addOption("dt_save", "Time step for saving.");
	OptionParser::addOption("tf", "Final time");
	OptionParser::addOption("dx", "Spatial discretization step in the x-direction (cm)");
	OptionParser::addOption("dy", "Spatial discretization step in the y-direction (cm)");
	OptionParser::addOption("Lx", "Horizontal length in cm");
	OptionParser::addOption("Ly", "Vertical length in cm");
	OptionParser::addOption("threads", "Number of threads");

	OptionParser::parseOptions(argc, argv);

	int model_index = OptionParser::parseInt("model");
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
	} else {
		return -1;
	}

	ODEMethod* method = new RushLarsenMethod();

	Nx = Lx / dx;
	Ny = Ly / dy;

	double* Y_old_ = new double [Nx * Ny * model->nStates]; 
	double* Y_new_ = new double [Nx * Ny * model->nStates];
	double* ALGS = new double [Nx * Ny * model->nAlgs]; double* PARAMS = new double [Nx * Ny * model->nParams];
	double* RHS = new double [Nx * Ny * (model->nStates + model->nStates_HH)];

	double t = 0, t_save = 0;

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

					// Stimulus
					double x = xi * dx, y = yi * dy;
					params[0] = (t > 1 && t < 3 && x > 0 && x < 0.1) ? 1 : -1;

					method->step(y_new_, model, params, algs, rhs, y_old_, t, dt);

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
					cout << "t = " << t << ", V(x0) = " << Y_old_[0] << ", V(xf) = " << Y_old_[(Nx * Ny - 1) * model->nStates] << endl;
					t_save = 0;
				}
			}
		}
	}

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;

    std::cout << "Simulation finished in " << elapsed.count() << " s." << std::endl;
}