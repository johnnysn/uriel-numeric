#include "pch.h"
#include "model/electrophysiology/Noble1962.h"
#include "model/electrophysiology/Noble1962_defs.h"

Noble1962::Noble1962() : CellModel(1, 3, 0, 0, 19, 10) {

}

void Noble1962::calc_algs_nl(double* algs, double* pars, double* Y_old_, double time)
{
	calc_i_Stim = calc_stimulus(pars, time);
    calc_g_Na = pow(Y_old_[1], 3.00000)*Y_old_[2]*g_Na_max;
    calc_g_K1 = 1.2*exp((((-Y_old_[0])-9.0e+01)/5.0e+01)) + (1.5e-02*exp(((Y_old_[0]+9.0e+01)/6.0e+01)));
    calc_g_K2 = 1.2*pow(Y_old_[3],4.0e+00);
    calc_i_Na = (calc_g_Na+1.4e-01)*(Y_old_[0] - E_Na);
	calc_i_K = (calc_g_K1+calc_g_K2)*(Y_old_[0]+100.000);
    calc_i_Leak = g_L*(Y_old_[0] - E_L);
}

void Noble1962::calc_algs_hh(double* algs, double* pars, double* Y_old_, double time)
{
	calc_alpha_m = (((1.0e-01*((-Y_old_[0])-4.8e+01))/(exp((((-Y_old_[0])-4.8e+01)/1.5e+01))-1.0e+00)));
	calc_beta_m = (((1.2e-01*(Y_old_[0]+8.0e+00))/(exp(((Y_old_[0]+8.0e+00)/5.0e+00))-1.0e+00)));
	calc_m_inf = calc_alpha_m / (calc_alpha_m + calc_beta_m);
	calc_tau_m = 1.000 / (calc_alpha_m + calc_beta_m);

	calc_alpha_h = ((1.7e-01*exp((((-Y_old_[0])-9.0e+01)/2.0e+01))));
    calc_beta_h = ((1.0/(1.0e+00+exp((((-Y_old_[0])-4.2e+01)/1.0e+01)))));
	calc_h_inf = calc_alpha_h / (calc_alpha_h + calc_beta_h);
	calc_tau_h = 1.000 / (calc_alpha_h + calc_beta_h);

	calc_alpha_n = (((1.0e-04*((-Y_old_[0])-5.0e+01))/(exp((((-Y_old_[0])-5.0e+01)/1.0e+01))-1.0e+00)));
    calc_beta_n = ((2.0e-03*exp((((-Y_old_[0])-9.0e+01)/8.0e+01))));
    calc_n_inf = calc_alpha_n / (calc_alpha_n + calc_beta_n);
	calc_tau_n = 1.000 / (calc_alpha_n + calc_beta_n);
}

void Noble1962::calc_rhs_nl(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
	calc_algs_nl(algs, pars, Y_old_, t);

	V_f_ = -(calc_i_Na + calc_i_K + calc_i_Leak)/Cm;
}

void Noble1962::calc_rhs_hh(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
	calc_algs_hh(algs, pars, Y_old_, t);

	m_f_ = (((calc_m_inf - m_old_) / calc_tau_m));
	h_f_ = (((calc_h_inf - h_old_) / calc_tau_h));
	n_f_ = (((calc_n_inf - n_old_) / calc_tau_n));
}

void Noble1962::calc_rhs_mk(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
}

void Noble1962::calc_hh_coeff(double* a, double* b, double* pars, double* algs, double* Y_old_, double t)
{
	calc_algs_hh(algs, pars, Y_old_, t);

	m_a_ = -1.0 / calc_tau_m;
	h_a_ = -1.0 / calc_tau_h;
	n_a_ = -1.0 / calc_tau_n;

	m_b_ = (((calc_m_inf) / calc_tau_m));
	h_b_ = (((calc_h_inf) / calc_tau_h));
	n_b_ = (((calc_n_inf) / calc_tau_n));
}

void Noble1962::set_default_parameters(double* pars)
{
	stim_state = 0;
	stim_amplitude = -5.20e+01;
	stim_period = 1.0e+03;
	stim_start = 0.0e+00;
	stim_duration = 0.0e+00;
	Cm = 12.0;
	g_Na_max = 400.0;
	E_Na = 40.0;
	g_L = 0.075;
	E_L = -60.0;
}

void Noble1962::set_default_initial_state(double* Y_old_)
{
	V_old_ = -87.0;
	m_old_ = 0.01;
	h_old_ = 0.8;
	n_old_ = 0.01;
}

double Noble1962::calc_stimulus(double* pars, double t)
{
	if (stim_state < 0)
		return 0;
	if (stim_state > 0)
		return stim_amplitude;

	double t_since_last_tick = t - floor(t / stim_period)*stim_period;
	double pulse_end = stim_start + stim_duration;
	if (t_since_last_tick >= stim_start && t_since_last_tick <= pulse_end) {
		return stim_amplitude;
	}
	else return 0;
}

void Noble1962::prep_mk_transitions(double* algs, double* pars, double* Y_old_, double t) {}
void Noble1962::calc_mk_transitions(double** Tr, int mk_index, double* pars, double* algs, double* Y_old_, double t) {}
