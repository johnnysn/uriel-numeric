#include "pch.h"
#include "model/electrophysiology/Fox2002.h"
#include "model/electrophysiology/Fox2002_defs.h"

Fox2002::Fox2002(): CellModel(3, 10, 0, 0, 52, 51)
{

}

void Fox2002::set_default_parameters(double* pars) 
{
	stim_state = 0;
	stim_amplitude = -8.0000000000e+01;
	stim_start = 5;
	stim_end = 10000.0e+03;
	stim_period = 400;
	stim_duration = 1.0000000000e+00;
	R = 8.3140000000e+00;
	T = 3.1000000000e+02;
	F = 9.6500000000e+01;
	Na_o = 1.3800000000e+02;
	Na_i = 1.0000000000e+01;
	g_Na = 1.2800000000e+01;
	shift_h = 0.0000000000e+00;
	shift_j = 0.0000000000e+00;
	g_K1 = 2.8000000000e+00;
	K_o = 4.0000000000e+00;
	K_mK1 = 1.3000000000e+01;
	K_i = 1.4940000000e+02;
	g_Kr = 1.3600000000e-02;
	g_Ks = 2.4500000000e-02;
	g_to = 2.3815000000e-01;
	g_Kp = 2.2160000000e-03;
	i_NaK_max = 6.9300000000e-01;
	K_mNai = 1.0000000000e+01;
	K_mKo = 1.5000000000e+00;
	K_NaCa = 1.5000000000e+03;
	K_mNa = 8.7500000000e+01;
	K_mCa = 1.3800000000e+03;
	Ca_o = 2.0000000000e+03;
	K_sat = 2.0000000000e-01;
	eta = 3.5000000000e-01;
	i_pCa_max = 5.0000000000e-02;
	K_mpCa = 5.0000000000e-02;
	g_Cab = 3.8420000000e-04;
	g_Nab = 3.1000000000e-03;
	P_Ca = 2.2600000000e-05;
	C_sc = 1.0000000000e+00;
	P_CaK = 5.7900000000e-07;
	i_Ca_half = -2.6500000000e-01;
	K_mfCa = 1.8000000000e-01;
	V_up = 1.0000000000e-01;
	K_mup = 3.2000000000e-01;
	P_rel = 6.0000000000e+00;
	P_leak = 1.0000000000e-06;
	CSQN_tot = 1.0000000000e+04;
	K_mCSQN = 6.0000000000e+02;
	V_myo = 2.5840000000e-05;
	V_SR = 2.0000000000e-06;
	CMDN_tot = 1.0000000000e+01;
	K_mCMDN = 2.0000000000e+00;
	A_Cap = 1.5340000000e-04;
}

void Fox2002::set_default_initial_state(double* Y_old_)
{
	V_old_ = -9.47000000e+01;
	m_old_ = 2.46760000e-04;
	h_old_ = 9.98690000e-01;
	j_old_ = 9.98870000e-01;
	X_kr_old_ = 2.29000000e-01;
	X_ks_old_ = 1.00000000e-04;
	X_to_old_ = 3.74200000e-05;
	Y_to_old_ = 1.00000000e+00;
	f_old_ = 9.83000000e-01;
	d_old_ = 1.00000000e-04;
	f_Ca_old_ = 9.42000000e-01;
	Ca_SR_old_ = 3.20000000e+02;
	Ca_i_old_ = 4.72000000e-02;
}

double Fox2002::calc_stimulus(double* pars, double t)
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

void Fox2002::calc_rhs_nl(double* Y_f_, double* pars, double* algs, double* Y_old_, double t) 
{
	calc_algs_nl(algs, pars, Y_old_, t);

	V_f_ = (-(calc_i_Na + calc_i_Ca + calc_i_CaK + calc_i_Kr + calc_i_Ks + calc_i_to + calc_i_K1 + calc_i_Kp + calc_i_NaCa + calc_i_NaK + calc_i_p_Ca + calc_i_Na_b + calc_i_Ca_b + calc_i_Stim));
	Ca_SR_f_ = ((calc_beta_SR*((calc_J_up - calc_J_leak) - calc_J_rel)*V_myo) / V_SR);
	Ca_i_f_ = (calc_beta_i*(((calc_J_rel + calc_J_leak) - calc_J_up) - (((A_Cap*C_sc) / (2.0e+00*F*V_myo))*((calc_i_Ca + calc_i_Ca_b + calc_i_p_Ca) - (2.0e+00*calc_i_NaCa)))));
}

void Fox2002::calc_rhs_hh(double* Y_f_, double* pars, double* algs, double* Y_old_, double t) 
{
	calc_algs_hh(algs, pars, Y_old_, t);

	m_f_ = ((calc_alpha_m*(1.0e+00 - m_old_)) - (calc_beta_m*m_old_));
	h_f_ = ((calc_alpha_h*(1.0e+00 - h_old_)) - (calc_beta_h*h_old_));
	j_f_ = ((calc_alpha_j*(1.0e+00 - j_old_)) - (calc_beta_j*j_old_));
	X_kr_f_ = ((calc_X_kr_inf - X_kr_old_) / calc_tau_X_kr);
	X_ks_f_ = ((calc_X_ks_infinity - X_ks_old_) / calc_tau_X_ks);
	X_to_f_ = ((calc_alpha_X_to*(1.0e+00 - X_to_old_)) - (calc_beta_X_to*X_to_old_));
	Y_to_f_ = ((calc_alpha_Y_to*(1.0e+00 - Y_to_old_)) - (calc_beta_Y_to*Y_to_old_));
	f_f_ = ((calc_f_infinity - f_old_) / calc_tau_f);
	d_f_ = ((calc_d_infinity - d_old_) / calc_tau_d);
	f_Ca_f_ = ((calc_f_Ca_infinity - f_Ca_old_) / calc_tau_f_Ca);
}

void Fox2002::calc_hh_coeff(double* as, double* bs, double* pars, double* algs, double* Y_old_, double t) 
{
	calc_algs_hh(algs, pars, Y_old_, t);

	m_a_ = -calc_alpha_m - calc_beta_m;
	h_a_ = -calc_alpha_h - calc_beta_h;
	j_a_ = -calc_alpha_j - calc_beta_j;
	X_kr_a_ = -1.0 / calc_tau_X_kr;
	X_ks_a_ = -1.0 / calc_tau_X_ks;
	X_to_a_ = -calc_alpha_X_to - calc_beta_X_to;
	Y_to_a_ = -calc_alpha_Y_to - calc_beta_Y_to;
	f_a_ = -1.0 / (calc_tau_f);
	d_a_ = -1.0 / calc_tau_d;
	f_Ca_a_ = -1.0 / calc_tau_f_Ca;

	m_b_ = calc_alpha_m;
	h_b_ = calc_alpha_h;
	j_b_ = calc_alpha_j;
	X_kr_b_ = calc_X_kr_inf / calc_tau_X_kr;
	X_ks_b_ = calc_X_ks_infinity / calc_tau_X_ks;
	X_to_b_ = calc_alpha_X_to;
	Y_to_b_ = calc_alpha_Y_to;
	f_b_ = calc_f_infinity / calc_tau_f;
	d_b_ = calc_d_infinity / calc_tau_d;
	f_Ca_b_ = calc_f_Ca_infinity / calc_tau_f_Ca;
}

void Fox2002::calc_algs_nl(double* algs, double* pars, double* Y_old_, double t) 
{
	calc_i_Stim = calc_stimulus(pars, t);
	calc_E_Na = (((R*T) / F)*log((Na_o / Na_i)));
	calc_E_K = (((R*T) / F)*log((K_o / K_i)));
	calc_R_V = (1.0e+00 / (1.0e+00 + (2.50e+00*exp((1.0e-01*(V_old_ + 2.80e+01))))));
	calc_E_Ks = (((R*T) / F)*log(((K_o + (1.8330e-02*Na_o)) / (K_i + (1.8330e-02*Na_i)))));
	calc_Kp_V = (1.0e+00 / (1.0e+00 + exp(((7.4880e+00 - V_old_) / 5.980e+00))));
	calc_sigma = ((1.0e+00 / 7.0e+00)*(exp((Na_o / 6.730e+01)) - 1.0e+00));
	calc_i_NaCa = ((K_NaCa / ((pow(K_mNa, 3.0e+00) + pow(Na_o, 3.0e+00))*(K_mCa + Ca_o)*(1.0e+00 + (K_sat*exp((((eta - 1.0e+00)*V_old_*F) / (R*T)))))))*((exp(((eta*V_old_*F) / (R*T)))*pow(Na_i, 3.0e+00)*Ca_o) - (exp((((eta - 1.0e+00)*V_old_*F) / (R*T)))*pow(Na_o, 3.0e+00)*Ca_i_old_)));
	calc_i_p_Ca = ((i_pCa_max*Ca_i_old_) / (K_mpCa + Ca_i_old_));
	calc_E_Ca = (((R*T) / (2.0e+00*F))*log((Ca_o / Ca_i_old_)));
	calc_i_Ca_max = (((((P_Ca / C_sc)*4.0e+00*V_old_*pow(F, 2.0e+00)) / (R*T))*((Ca_i_old_*exp(((2.0e+00*V_old_*F) / (R*T)))) - (3.410e-01*Ca_o))) / (exp(((2.0e+00*V_old_*F) / (R*T))) - 1.0e+00));
	calc_J_up = (V_up / (1.0e+00 + pow((K_mup / Ca_i_old_), 2.0e+00)));
	calc_gamma = (1.0e+00 / (1.0e+00 + pow((2.0e+03 / Ca_SR_old_), 3.0e+00)));
	calc_J_leak = (P_leak*(Ca_SR_old_ - Ca_i_old_));
	calc_beta_SR = (1.0e+00 / (1.0e+00 + ((CSQN_tot*K_mCSQN) / pow((K_mCSQN + Ca_SR_old_), 2.0e+00))));
	calc_beta_i = (1.0e+00 / (1.0e+00 + ((CMDN_tot*K_mCMDN) / pow((K_mCMDN + Ca_i_old_), 2.0e+00))));
	calc_i_Na = (g_Na*pow(m_old_, 3.0e+00)*h_old_*j_old_*(V_old_ - calc_E_Na));
	calc_i_Kr = (g_Kr*calc_R_V*X_kr_old_*pow((K_o / 4.0e+00), 1.0 / 2.0)*(V_old_ - calc_E_K));
	calc_i_Ks = (g_Ks*pow(X_ks_old_, 2.0e+00)*(V_old_ - calc_E_Ks));
	calc_i_to = (g_to*X_to_old_*Y_to_old_*(V_old_ - calc_E_K));
	calc_i_Ca_b = (g_Cab*(V_old_ - calc_E_Ca));
	calc_i_Na_b = (g_Nab*(V_old_ - calc_E_Na));
	calc_i_CaK = (((((((P_CaK / C_sc)*f_old_*d_old_*f_Ca_old_) / (1.0e+00 + (calc_i_Ca_max / i_Ca_half)))*1.0e+03*V_old_*pow(F, 2.0e+00)) / (R*T))*((K_i*exp(((V_old_*F) / (R*T)))) - K_o)) / (exp(((V_old_*F) / (R*T))) - 1.0e+00));
	calc_J_rel = ((P_rel*f_old_*d_old_*f_Ca_old_*((calc_gamma*Ca_SR_old_) - Ca_i_old_)) / (1.0e+00 + (1.650e+00*exp((V_old_ / 2.0e+01)))));
	calc_K1_infinity = (1.0e+00 / (2.0e+00 + exp((((1.620e+00*F) / (R*T))*(V_old_ - calc_E_K)))));
	calc_i_Kp = (g_Kp*calc_Kp_V*(V_old_ - calc_E_K));
	calc_f_NaK = (1.0e+00 / (1.0e+00 + (1.2450e-01*exp((((-1.0e-01)*V_old_*F) / (R*T)))) + (3.650e-02*calc_sigma*exp((((-V_old_)*F) / (R*T))))));
	calc_i_Ca = (calc_i_Ca_max*f_old_*d_old_*f_Ca_old_);
	calc_i_NaK = ((((i_NaK_max*calc_f_NaK) / (1.0e+00 + pow((K_mNai / Na_i), 1.50e+00)))*K_o) / (K_o + K_mKo));
	calc_i_K1 = (((g_K1*calc_K1_infinity*K_o) / (K_o + K_mK1))*(V_old_ - calc_E_K));
}

void Fox2002::calc_algs_hh(double* algs, double* pars, double* Y_old_, double time) 
{
	calc_E0_m = (V_old_ + 4.7130e+01);
	calc_beta_m = (8.0e-02*exp(((-V_old_) / 1.10e+01)));
	calc_alpha_m = ((3.20e-01*calc_E0_m) / (1.0e+00 - exp(((-1.0e-01)*calc_E0_m))));
	calc_alpha_h = (1.350e-01*exp((((V_old_ + 8.0e+01) - shift_h) / (-6.80e+00))));
	calc_beta_h = (7.50e+00 / (1.0e+00 + exp(((-1.0e-01)*((V_old_ + 1.10e+01) - shift_h)))));
	calc_alpha_j = ((1.750e-01*exp((((V_old_ + 1.0e+02) - shift_j) / (-2.30e+01)))) / (1.0e+00 + exp((1.50e-01*((V_old_ + 7.90e+01) - shift_j)))));
	calc_beta_j = (3.0e-01 / (1.0e+00 + exp(((-1.0e-01)*((V_old_ + 3.20e+01) - shift_j)))));
	calc_X_kr_inf = (1.0e+00 / (1.0e+00 + exp(((-2.1820e+00) - (1.8190e-01*V_old_)))));
	calc_tau_X_kr = (4.30e+01 + (1.0e+00 / (exp(((-5.4950e+00) + (1.6910e-01*V_old_))) + exp(((-7.6770e+00) - (1.280e-02*V_old_))))));
	calc_X_ks_infinity = (1.0e+00 / (1.0e+00 + exp(((V_old_ - 1.60e+01) / (-1.360e+01)))));
	calc_tau_X_ks = (1.0e+00 / (((7.190e-05*(V_old_ - 1.0e+01)) / (1.0e+00 - exp(((-1.480e-01)*(V_old_ - 1.0e+01))))) + ((1.310e-04*(V_old_ - 1.0e+01)) / (exp((6.870e-02*(V_old_ - 1.0e+01))) - 1.0e+00))));
	calc_alpha_X_to = (4.5160e-02*exp((3.5770e-02*V_old_)));
	calc_beta_X_to = (9.890e-02*exp(((-6.2370e-02)*V_old_)));
	calc_alpha_Y_to = ((5.4150e-03*exp(((V_old_ + 3.350e+01) / (-5.0e+00)))) / (1.0e+00 + (5.13350e-02*exp(((V_old_ + 3.350e+01) / (-5.0e+00))))));
	calc_beta_Y_to = ((5.4150e-03*exp(((V_old_ + 3.350e+01) / 5.0e+00))) / (1.0e+00 + (5.13350e-02*exp(((V_old_ + 3.350e+01) / 5.0e+00)))));
	calc_tau_f = (3.0e+01 + 2.0e+02 / (1.0e+00 + exp((V_old_ + 2.0e+01) / 9.50e+00)));
	calc_f_infinity = (1.0e+00 / (1.0e+00 + exp(((V_old_ + 1.250e+01) / 5.0e+00))));
	calc_d_infinity = (1.0e+00 / (1.0e+00 + exp(((V_old_ + 1.0e+01) / (-6.240e+00)))));
	calc_E0_m_duplicated_L_type_Ca_current_d_gate = (V_old_ + 4.0e+01);
	calc_tau_d = (1.0e+00 / (((2.50e-01*exp(((-1.0e-02)*V_old_))) / (1.0e+00 + exp(((-7.0e-02)*V_old_)))) + ((7.0e-02*exp(((-5.0e-02)*calc_E0_m_duplicated_L_type_Ca_current_d_gate))) / (1.0e+00 + exp((5.0e-02*calc_E0_m_duplicated_L_type_Ca_current_d_gate))))));
	calc_f_Ca_infinity = (1.0e+00 / (1.0e+00 + pow((Ca_i_old_ / K_mfCa), 3.0e+00)));
	calc_tau_f_Ca = 3.0e+01;
}