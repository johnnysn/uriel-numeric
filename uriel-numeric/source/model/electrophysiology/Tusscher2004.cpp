#include "pch.h"
#include "model/electrophysiology/Tusscher2004.h"
#include "model/electrophysiology/Tusscher2004_defs.h"

Tusscher2004::Tusscher2004() : CellModel(5, 12, 0, 0, 68, 46) {

}

void Tusscher2004::calc_algs_nl(double* algs, double* pars, double* Y_old_, double time)
{
	calc_i_Stim = calc_stimulus(pars, time);	//0
	calc_E_Na = (((R*T) / F)*log((Na_o / Na_i_old_)));	//2
	calc_E_K = (((R*T) / F)*log((K_o / K_i_old_)));	//3
	calc_E_Ks = (((R*T) / F)*log(((K_o + (P_kna*Na_o)) / (K_i_old_ + (P_kna*Na_i_old_)))));	//4
	calc_E_Ca = (((5.0e-01*R*T) / F)*log((Ca_o / Ca_i_old_)));	//5
	calc_i_CaL = ((((g_CaL*d_old_*f_old_*fCa_old_*4.0e+00*V_old_*pow(F, 2.0e+00)) / (R*T))*((Ca_i_old_*exp(((2.0e+00*V_old_*F) / (R*T)))) - (3.410e-01*Ca_o))) / (exp(((2.0e+00*V_old_*F) / (R*T))) - 1.0e+00));	//44
	calc_i_NaK = (((((P_NaK*K_o) / (K_o + K_mk))*Na_i_old_) / (Na_i_old_ + K_mNa)) / (1.0e+00 + (1.2450e-01*exp((((-1.0e-01)*V_old_*F) / (R*T)))) + (3.530e-02*exp((((-V_old_)*F) / (R*T))))));	//69
	calc_i_NaCa = ((K_NaCa*((exp(((gamma*V_old_*F) / (R*T)))*pow(Na_i_old_, 3.0e+00)*Ca_o) - (exp((((gamma - 1.0e+00)*V_old_*F) / (R*T)))*pow(Na_o, 3.0e+00)*Ca_i_old_*alpha))) / ((pow(Km_Nai, 3.0e+00) + pow(Na_o, 3.0e+00))*(Km_Ca + Ca_o)*(1.0e+00 + (K_sat*exp((((gamma - 1.0e+00)*V_old_*F) / (R*T)))))));	//70
	calc_i_p_Ca = ((g_pCa*Ca_i_old_) / (Ca_i_old_ + K_pCa));	//71
	calc_i_rel = ((((a_rel*pow(Ca_SR_old_, 2.0e+00)) / (pow(b_rel, 2.0e+00) + pow(Ca_SR_old_, 2.0e+00))) + c_rel)*d_old_*g_old_);	//73
	calc_i_up = (Vmax_up / (1.0e+00 + (pow(K_up, 2.0e+00) / pow(Ca_i_old_, 2.0e+00))));	//74
	calc_i_leak = (V_leak*(Ca_SR_old_ - Ca_i_old_));	//75
	calc_Ca_i_bufc = (1.0e+00 / (1.0e+00 + ((Buf_c*K_buf_c) / pow((Ca_i_old_ + K_buf_c), 2.0e+00))));	//79
	calc_Ca_sr_bufsr = (1.0e+00 / (1.0e+00 + ((Buf_sr*K_buf_sr) / pow((Ca_SR_old_ + K_buf_sr), 2.0e+00))));	//80
	calc_i_Kr = (g_Kr*pow((K_o / 5.40e+00), 1.0 / 2.0)*Xr1_old_*Xr2_old_*(V_old_ - calc_E_K));	//10
	calc_i_Ks = (g_Ks*pow(Xs_old_, 2.0e+00)*(V_old_ - calc_E_Ks));	//21
	calc_i_Na = (g_Na*pow(m_old_, 3.0e+00)*h_old_*j_old_*(V_old_ - calc_E_Na));	//27
	calc_i_b_Na = (g_bna*(V_old_ - calc_E_Na));	//43
	calc_i_b_Ca = (g_bca*(V_old_ - calc_E_Ca));	//61
	calc_i_to = (g_to*r_old_*s_old_*(V_old_ - calc_E_K));	//62
	calc_i_p_K = ((g_pK*(V_old_ - calc_E_K)) / (1.0e+00 + exp(((2.50e+01 - V_old_) / 5.980e+00))));	//72
	calc_alpha_K1 = (1.0e-01 / (1.0e+00 + exp((6.0e-02*((V_old_ - calc_E_K) - 2.0e+02)))));	//6
	calc_beta_K1 = (((3.0e+00*exp((2.0e-04*((V_old_ - calc_E_K) + 1.0e+02)))) + (1.0e+00*exp((1.0e-01*((V_old_ - calc_E_K) - 1.0e+01))))) / (1.0e+00 + exp(((-5.0e-01)*(V_old_ - calc_E_K)))));	//7
	calc_xK1_inf = (calc_alpha_K1 / (calc_alpha_K1 + calc_beta_K1));	//8
	calc_i_K1 = (g_K1*calc_xK1_inf*pow((K_o / 5.40e+00), 1.0 / 2.0)*(V_old_ - calc_E_K));	//9
}

void Tusscher2004::calc_algs_hh(double* algs, double* pars, double* Y_old_, double time)
{
	calc_xr1_inf = (1.0e+00 / (1.0e+00 + exp((((-2.60e+01) - V_old_) / 7.0e+00))));	//11
	calc_alpha_xr1 = (4.50e+02 / (1.0e+00 + exp((((-4.50e+01) - V_old_) / 1.0e+01))));	//12
	calc_beta_xr1 = (6.0e+00 / (1.0e+00 + exp(((V_old_ + 3.0e+01) / 1.150e+01))));	//13
	calc_xr2_inf = (1.0e+00 / (1.0e+00 + exp(((V_old_ + 8.80e+01) / 2.40e+01))));	//16
	calc_alpha_xr2 = (3.0e+00 / (1.0e+00 + exp((((-6.0e+01) - V_old_) / 2.0e+01))));	//17
	calc_beta_xr2 = (1.120e+00 / (1.0e+00 + exp(((V_old_ - 6.0e+01) / 2.0e+01))));	//18
	calc_xs_inf = (1.0e+00 / (1.0e+00 + exp((((-5.0e+00) - V_old_) / 1.40e+01))));	//22
	calc_alpha_xs = (1.10e+03 / pow((1.0e+00 + exp((((-1.0e+01) - V_old_) / 6.0e+00))), 1.0 / 2.0));	//23
	calc_beta_xs = (1.0e+00 / (1.0e+00 + exp(((V_old_ - 6.0e+01) / 2.0e+01))));	//24

	calc_m_inf = (1.0e+00 / pow((1.0e+00 + exp((((-5.6860e+01) - V_old_) / 9.03e+00))), 2.0e+00));	//28
	calc_alpha_m = (1.0e+00 / (1.0e+00 + exp((((-6.0e+01) - V_old_) / 5.0e+00))));	//29
	calc_beta_m = ((1.0e-01 / (1.0e+00 + exp(((V_old_ + 3.50e+01) / 5.0e+00)))) + (1.0e-01 / (1.0e+00 + exp(((V_old_ - 5.0e+01) / 2.0e+02)))));	//30
	calc_h_inf = (1.0e+00 / pow((1.0e+00 + exp(((V_old_ + 7.1550e+01) / 7.430e+00))), 2.0e+00));	//33
	calc_alpha_h = (V_old_ < -40) ? 5.70e-02*exp(((-(V_old_ + 8.0e+01)) / 6.80e+00)) : 0;	//34
	calc_beta_h = (V_old_ < -40) ? (2.70e+00*exp((7.90e-02*V_old_))) + (3.10e+05*exp((3.4850e-01*V_old_))) : 7.70e-01 / (1.30e-01*(1.0 + exp(((V_old_ + 1.0660e+01) / (-1.110e+01)))));	//35
	calc_j_inf = (1.0e+00 / pow((1.0e+00 + exp(((V_old_ + 7.1550e+01) / 7.430e+00))), 2.0e+00));	//38
	calc_alpha_j = (V_old_ < -40) ? ((-2.5428e4)*exp(0.2444*V_old_) - (6.948e-6)*exp(-0.04391*V_old_))*(V_old_ + 37.78) / (1. + exp(0.311*(V_old_ + 79.23))) : 0;	//39
	calc_beta_j = (V_old_ < -40) ? 0.02424*exp(-0.01052*V_old_) / (1. + exp(-0.1378*(V_old_ + 40.14))) : 0.6*exp((0.057)*V_old_) / (1. + exp(-0.1*(V_old_ + 32.)));	//40

	calc_d_inf = (1.0e+00 / (1.0e+00 + exp((((-5.0e+00) - V_old_) / 7.50e+00))));	//45
	calc_alpha_d = ((1.40e+00 / (1.0e+00 + exp((((-3.50e+01) - V_old_) / 1.30e+01)))) + 2.50e-01);	//46
	calc_beta_d = (1.40e+00 / (1.0e+00 + exp(((V_old_ + 5.0e+00) / 5.0e+00))));	//47
	calc_gamma_d = (1.0e+00 / (1.0e+00 + exp(((5.0e+01 - V_old_) / 2.0e+01))));	//48

	calc_f_inf = 1. / (1. + exp((V_old_ + 20) / 7));	//51
	calc_tau_f = 1125 * exp(-(V_old_ + 27)*(V_old_ + 27) / 240) + 80 + 165 / (1. + exp((25 - V_old_) / 10));	//52 300 -> 240 ?
	//calc_tau_f = ( 3.0e+01+3.5e+02/(1.0e+00+exp((V_on_f+2.5e+01)/9.50e+00)) );

	calc_tau_fCa = 2.0;	//58
	calc_s_inf = (1.0e+00 / (1.0e+00 + exp(((V_old_ + 2.0e+01) / 5.0e+00))));	//63
	calc_tau_s = ((8.50e+01*exp(((-pow((V_old_ + 4.50e+01), 2.0e+00)) / 3.20e+02))) + (5.0e+00 / (1.0e+00 + exp(((V_old_ - 2.0e+01) / 5.0e+00)))) + 3.0e+00);	//64
	calc_r_inf = (1.0e+00 / (1.0e+00 + exp(((2.0e+01 - V_old_) / 6.0e+00))));	//66
	calc_tau_r = ((9.50e+00*exp(((-pow((V_old_ + 4.0e+01), 2.0e+00)) / 1.80e+03))) + 8.0e-01);	//67
	calc_g_inf = Ca_i_old_ < 3.50e-04 ? 1.0e+00 / (1.0e+00 + pow((Ca_i_old_ / 3.50e-04), 6.0e+00)) : 1.0e+00 / (1.0e+00 + pow((Ca_i_old_ / 3.50e-04), 1.60e+01));	//76
	calc_tau_xr1 = (1.0e+00*calc_alpha_xr1*calc_beta_xr1);	//14
	calc_tau_xr2 = (1.0e+00*calc_alpha_xr2*calc_beta_xr2);	//19
	calc_tau_xs = (1.0e+00*calc_alpha_xs*calc_beta_xs);	//25
	calc_tau_m = (1.0e+00*calc_alpha_m*calc_beta_m);	//31
	calc_tau_h = (1.0e+00 / (calc_alpha_h + calc_beta_h));	//36
	calc_tau_j = (1.0e+00 / (calc_alpha_j + calc_beta_j));	//41
	calc_tau_d = ((1.0e+00*calc_alpha_d*calc_beta_d) + calc_gamma_d);	//49

	calc_alpha_fCa = (1.0e+00 / (1.0e+00 + pow((Ca_i_old_ / 3.250e-04), 8.0e+00)));	//54
	calc_beta_fCa = (1.0e-01 / (1.0e+00 + exp(((Ca_i_old_ - 5.0e-04) / 1.0e-04))));	//55
	calc_gama_fCa = (2.0e-01 / (1.0e+00 + exp(((Ca_i_old_ - 7.50e-04) / 8.0e-04))));	//56
	calc_fCa_inf = ((calc_alpha_fCa + calc_beta_fCa + calc_gama_fCa + 2.30e-01) / 1.460e+00);	//57

	calc_d_g = ((calc_g_inf - g_old_) / tau_g);	//77
	calc_d_fCa = ((calc_fCa_inf - fCa_old_) / calc_tau_fCa);	//59
}

void Tusscher2004::calc_rhs_nl(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
	calc_algs_nl(algs, pars, Y_old_, t);

	V_f_ = -(calc_i_K1 + calc_i_to + calc_i_Kr + calc_i_Ks + calc_i_CaL + calc_i_NaK + calc_i_Na + calc_i_b_Na + calc_i_NaCa + calc_i_b_Ca + calc_i_p_K + calc_i_p_Ca + calc_i_Stim);
	Ca_i_f_ = ((calc_Ca_i_bufc*(((calc_i_leak - calc_i_up) + calc_i_rel) - (((1.0e+00*((calc_i_CaL + calc_i_b_Ca + calc_i_p_Ca) - (2.0e+00*calc_i_NaCa))) / (2.0e+00*1.0e+00*V_c*F))*Cm))));	// 81
	Ca_SR_f_ = ((((calc_Ca_sr_bufsr*V_c) / V_sr)*(calc_i_up - (calc_i_rel + calc_i_leak))));	// 82
	Na_i_f_ = ((((-1.0e+00)*(calc_i_Na + calc_i_b_Na + (3.0e+00*calc_i_NaK) + (3.0e+00*calc_i_NaCa))*Cm) / (1.0e+00*V_c*F)));	// 83
	K_i_f_ = ((((-1.0e+00)*((calc_i_K1 + calc_i_to + calc_i_Kr + calc_i_Ks + calc_i_p_K + calc_i_Stim) - (2.0e+00*calc_i_NaK))*Cm) / (1.0e+00*V_c*F)));	// 84
}

void Tusscher2004::calc_rhs_hh(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
	calc_algs_hh(algs, pars, Y_old_, t);

	Xr1_f_ = (((calc_xr1_inf - Xr1_old_) / calc_tau_xr1));	// 15
	Xr2_f_ = (((calc_xr2_inf - Xr2_old_) / calc_tau_xr2));	// 20
	Xs_f_ = (((calc_xs_inf - Xs_old_) / calc_tau_xs));	// 26
	m_f_ = (((calc_m_inf - m_old_) / calc_tau_m));	// 32
	h_f_ = (((calc_h_inf - h_old_) / calc_tau_h));	// 37
	j_f_ = (((calc_j_inf - j_old_) / calc_tau_j));	// 42
	d_f_ = (((calc_d_inf - d_old_) / calc_tau_d));	// 50
	f_f_ = (((calc_f_inf - f_old_) / calc_tau_f));	// 53
	fCa_f_ = (calc_fCa_inf > fCa_old_&& V_old_ > -6.0e+01) ? 0.0 : calc_d_fCa;	// 60
	s_f_ = (((calc_s_inf - s_old_) / calc_tau_s));	// 65
	r_f_ = (((calc_r_inf - r_old_) / calc_tau_r));	// 68
	g_f_ = (calc_g_inf > g_old_&&V_old_ > -6.0e+01) ? 0.0 : calc_d_g;	// 78
}

void Tusscher2004::calc_rhs_mk(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
}

void Tusscher2004::calc_hh_coeff(double* a, double* b, double* pars, double* algs, double* Y_old_, double t)
{
	calc_algs_hh(algs, pars, Y_old_, t);

	Xr1_a_ = -1.0 / calc_tau_xr1;	// 15
	Xr2_a_ = -1.0 / calc_tau_xr2;	// 20
	Xs_a_ = -1.0 / calc_tau_xs;	// 26
	m_a_ = -1.0 / calc_tau_m;	// 32
	h_a_ = -1.0 / calc_tau_h;	// 37
	j_a_ = -1.0 / calc_tau_j;	// 42
	d_a_ = -1.0 / calc_tau_d;	// 50
	f_a_ = -1.0 / calc_tau_f;	// 53
	//fCa_a_= (calc_fCa_inf>fCa_old_&& V_old_>-6.0e+01) ? 0.0 : -1.0/calc_tau_fCa;	// 56
	fCa_a_ = -1.0 / calc_tau_fCa;
	s_a_ = -1.0 / calc_tau_s;	// 64
	r_a_ = -1.0 / calc_tau_r;	// 67
	//g_a_=(calc_g_inf>g_old_&&V_old_>-6.0e+01) ? 0.0 : -1.0/tau_g;
	g_a_ = -1.0 / tau_g;

	Xr1_b_ = (((calc_xr1_inf) / calc_tau_xr1));
	Xr2_b_ = (((calc_xr2_inf) / calc_tau_xr2));
	Xs_b_ = (((calc_xs_inf) / calc_tau_xs));
	m_b_ = (((calc_m_inf) / calc_tau_m));
	h_b_ = (((calc_h_inf) / calc_tau_h));
	j_b_ = (((calc_j_inf) / calc_tau_j));
	d_b_ = (((calc_d_inf) / calc_tau_d));
	f_b_ = (((calc_f_inf) / calc_tau_f));
	//fCa_b_= (calc_fCa_inf>fCa_old_&& V_old_>-6.0e+01) ? 0.0 : calc_fCa_inf/calc_tau_fCa;
	fCa_b_ = calc_fCa_inf / calc_tau_fCa;
	s_b_ = (((calc_s_inf) / calc_tau_s));
	r_b_ = (((calc_r_inf) / calc_tau_r));
	//g_b_ = (calc_g_inf>g_old_&&V_old_>-6.0e+01) ? 0.0 : calc_g_inf/tau_g;
	g_b_ = calc_g_inf / tau_g;

	if ((fCa_old_ * fCa_a_ + fCa_b_) > 0.0 && V_old_ > -37) {
		fCa_a_ = fCa_b_ = 0.0;	// 56
	}
	if ((g_old_ * g_a_ + g_b_) > 0.0 && V_old_ > -37) {
		g_a_ = g_b_ = 0.0;	// 56
	}
}

void Tusscher2004::set_default_parameters(double* pars)
{
	stim_state = 0;
	stim_amplitude = -5.20e+01;
	stim_period = 1.0e+03;
	stim_start = 5.0e+00;
	stim_duration = 1.0e+00;
	R = 8.3144720e+03;
	T = 3.10e+02;
	F = 9.64853415e+04;
	Na_o = 1.40e+02;
	K_o = 5.40e+00;
	P_kna = 3.0e-02;
	Ca_o = 2.0e+00;
	g_K1 = 5.4050e+00;
	g_Kr = 0.096; //0.134
	g_Ks = 0.245; //0.270
	g_Na = 1.48380e+01;
	g_bna = 2.90e-04;
	g_CaL = 1.750e-04;
	g_bca = 5.920e-04;
	g_to = 2.940e-01;
	P_NaK = 1.3620e+00;
	K_mk = 1.0e+00;
	K_mNa = 4.0e+01;
	K_NaCa = 1.0e+03;
	gamma = 3.50e-01;
	alpha = 2.50e+00;
	Km_Nai = 8.750e+01;
	Km_Ca = 1.380e+00;
	K_sat = 1.0e-01;
	g_pCa = 8.250e-01;
	K_pCa = 5.0e-04;
	g_pK = 1.460e-02;
	a_rel = 1.64640e-02;
	b_rel = 2.50e-01;
	c_rel = 8.2320e-03;
	Vmax_up = 4.250e-04;
	K_up = 2.50e-04;
	V_leak = 8.0e-05;
	tau_g = 2.0e+00;
	Buf_c = 1.50e-01;
	K_buf_c = 1.0e-03;
	Buf_sr = 1.0e+01;
	K_buf_sr = 3.0e-01;
	V_c = 1.64040e-02;
	Cm = 1.850e-01;
	V_sr = 1.0940e-03;
}

void Tusscher2004::set_default_initial_state(double* Y_old_)
{
	V_old_ = -8.620e+01;
	Xr1_old_ = 0.0e+00;
	Xr2_old_ = 1.0e+00;
	Xs_old_ = 0.0e+00;
	m_old_ = 0.0e+00;
	h_old_ = 7.50e-01;
	j_old_ = 7.50e-01;
	d_old_ = 0.0e+00;
	f_old_ = 1.0e+00;
	fCa_old_ = 1.0e+00;
	s_old_ = 1.0e+00;
	r_old_ = 0.0e+00;
	g_old_ = 1.0e+00;
	Ca_i_old_ = 0.00008;//2.0e-04;
	Ca_SR_old_ = 0.56;//2.0e-01;
	Na_i_old_ = 11.6;//1.160e+01;
	K_i_old_ = 138.3;//1.3830e+02;
}

double Tusscher2004::calc_stimulus(double* pars, double t)
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