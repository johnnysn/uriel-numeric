#include "pch.h"
#include "model/electrophysiology/Bondarenko2004.h"
#include "model/electrophysiology/Bondarenko2004_defs.h"

Bondarenko2004::Bondarenko2004() : CellModel(11, 8, 26, 4, 68, 72)
{
    nStates_MKM[0] = 9;
    nStates_MKM[1] = 8;
    nStates_MKM[2] = 5;
    nStates_MKM[3] = 4;

    nStates_MKM_max = 9;
}

void Bondarenko2004::set_default_parameters(double* pars)
{
    stim_state = 0;
    stim_amplitude = -8.0e+01;
    stim_start = 1.0e+0;
    stim_end = 1.0e+05;
    stim_period = 60.0;
    stim_duration = 5.0e-01;
    Acap = 1.5340e-04;
    Cm = 1.0e+00;
    Vmyo = 2.5840e-05;
    F = 9.650e+01;
    VJSR = 1.20e-07;
    Vss = 1.4850e-09;
    VNSR = 2.0980e-06;
    CMDN_tot = 5.0e+01;
    Km_CMDN = 2.380e-01;
    CSQN_tot = 1.50e+04;
    Km_CSQN = 8.0e+02;
    v1 = 4.50e+00;
    tau_tr = 2.0e+01;
    tau_xfer = 8.0e+00;
    v2 = 1.740e-05;
    v3 = 4.50e-01;
    Km_up = 5.0e-01;
    k_plus_htrpn = 2.370e-03;
    HTRPN_tot = 1.40e+02;
    k_plus_ltrpn = 3.270e-02;
    LTRPN_tot = 7.0e+01;
    k_minus_htrpn = 3.20e-05;
    k_minus_ltrpn = 1.960e-02;
    i_CaL_max = 7.0e+00;
    k_plus_a = 6.0750e-03;
    n = 4.0e+00;
    k_minus_b = 9.650e-01;
    k_minus_c = 8.0e-04;
    k_minus_a = 7.1250e-02;
    k_plus_b = 4.050e-03;
    m = 3.0e+00;
    k_plus_c = 9.0e-03;
    g_CaL = 1.7290e-01;
    E_CaL = 6.30e+01;
    Kpcb = 5.0e-04;
    Kpc_max = 2.33240e-01;
    Kpc_half = 2.0e+01;
    i_pCa_max = 1.0e+00;
    Km_pCa = 5.0e-01;
    k_NaCa = 2.9280e+02;
    K_mNa = 8.750e+04;
    Nao = 1.40e+05;
    K_mCa = 1.380e+03;
    Cao = 1.80e+03;
    k_sat = 1.0e-01;
    eta = 3.50e-01;
    R = 8.3140e+00;
    T = 2.980e+02;
    g_Cab = 3.670e-04;
    g_Na = 1.30e+01;
    Ko = 5.40e+03;
    g_Nab = 2.60e-03;
    g_Kto_f = 4.0670e-01;
    g_Kto_s = 0.0e+00;
    g_Ks = 5.750e-03;
    g_Kur = 1.60e-01;
    g_Kss = 5.0e-02;
    g_Kr = 7.80e-02;
    kf = 2.37610e-02;
    kb = 3.67780e-02;
    i_NaK_max = 8.80e-01;
    Km_Nai = 2.10e+04;
    Km_Ko = 1.50e+03;
    g_ClCa = 1.0e+01;
    Km_Cl = 1.0e+01;
    E_Cl = -4.0e+01;
}

void Bondarenko2004::set_default_initial_state(double* Y_old_)
{
    V_old_ = -8.242020e+01;
    Cai_old_ = 1.150010e-01;
    Cass_old_ = 1.150010e-01;
    CaJSR_old_ = 1.29950e+03;
    CaNSR_old_ = 1.29950e+03;
    P_RyR_old_ = 0.0e+00;
    LTRPN_Ca_old_ = 1.126840e+01;
    HTRPN_Ca_old_ = 1.25290e+02;
    P_O1_old_ = 1.491020e-05;
    P_O2_old_ = 9.517260e-11;
    P_C2_old_ = 1.67740e-04;
    O_old_ = 9.303080e-19;
    C2_old_ = 1.242160e-04;
    C3_old_ = 5.786790e-09;
    C4_old_ = 1.198160e-13;
    I1_old_ = 4.979230e-19;
    I2_old_ = 3.458470e-14;
    I3_old_ = 1.851060e-14;
    Nai_old_ = 1.423710e+04;
    C_Na2_old_ = 2.07520e-02;
    C_Na1_old_ = 2.791320e-04;
    O_Na_old_ = 7.134830e-07;
    IF_Na_old_ = 1.531760e-04;
    I1_Na_old_ = 6.733450e-07;
    I2_Na_old_ = 1.557870e-09;
    IC_Na2_old_ = 1.138790e-02;
    IC_Na3_old_ = 3.42780e-01;
    Ki_old_ = 1.43720e+05;
    ato_f_old_ = 2.655630e-03;
    ito_f_old_ = 9.999770e-01;
    ato_s_old_ = 4.170690e-04;
    ito_s_old_ = 9.985430e-01;
    nKs_old_ = 2.627530e-04;
    aur_old_ = 4.170690e-04;
    iur_old_ = 9.985430e-01;
    aKss_old_ = 4.170690e-04;
    iKss_old_ = 1.0e+00;
    C_K2_old_ = 6.412290e-04;
    C_K1_old_ = 9.925130e-04;
    O_K_old_ = 1.752980e-04;
    I_K_old_ = 3.191290e-05;

    P_C1_old_ = (1.0e+00 - (P_C2_old_ + P_O1_old_ + P_O2_old_));	//19
    C1_old_ = (1.0e+00 - (O_old_ + C2_old_ + C3_old_ + C4_old_ + I1_old_ + I2_old_ + I3_old_));	//24
    C_Na3_old_ = (1.0e+00 - (O_Na_old_ + C_Na1_old_ + C_Na2_old_ + IF_Na_old_ + I1_Na_old_ + I2_Na_old_ + IC_Na2_old_ + IC_Na3_old_));
    C_K0_old_ = (1.0e+00 - (C_K1_old_ + C_K2_old_ + O_K_old_ + I_K_old_));
}

double Bondarenko2004::calc_stimulus(double* pars, double t) 
{
    if (stim_state < 0)
        return 0;
    if (stim_state > 0)
        return stim_amplitude;

    double t_since_last_tick = t - floor(t / stim_period) * stim_period;
    double pulse_end = stim_start + stim_duration;
    if (t_since_last_tick >= stim_start && t_since_last_tick <= pulse_end) {
        return stim_amplitude;
    }
    else return 0;
}

void Bondarenko2004::calc_algs_nl(double* algs, double* pars, double* Y_old_, double time)
{
    calc_i_stim = calc_stimulus(pars, time);	//0
    calc_Bi = pow((1.0e+00 + ((CMDN_tot * Km_CMDN) / pow((Km_CMDN + Cai_old_), 2.0e+00))), (-1.0e+00));	//6
    calc_Bss = pow((1.0e+00 + ((CMDN_tot * Km_CMDN) / pow((Km_CMDN + Cass_old_), 2.0e+00))), (-1.0e+00));	//7
    calc_BJSR = pow((1.0e+00 + ((CSQN_tot * Km_CSQN) / pow((Km_CSQN + CaJSR_old_), 2.0e+00))), (-1.0e+00));	//8
    calc_J_rel = (v1 * (P_O1_old_ + P_O2_old_) * (CaJSR_old_ - Cass_old_) * P_RyR_old_);	//9
    calc_J_tr = ((CaNSR_old_ - CaJSR_old_) / tau_tr);	//10
    calc_J_xfer = ((Cass_old_ - Cai_old_) / tau_xfer);	//11
    calc_J_leak = (v2 * (CaNSR_old_ - Cai_old_));	//12
    calc_J_up = ((v3 * pow(Cai_old_, 2.0e+00)) / (pow(Km_up, 2.0e+00) + pow(Cai_old_, 2.0e+00)));	//13
    calc_J_trpn = (((k_plus_htrpn * Cai_old_ * (HTRPN_tot - HTRPN_Ca_old_)) + (k_plus_ltrpn * Cai_old_ * (LTRPN_tot - LTRPN_Ca_old_))) - ((k_minus_htrpn * HTRPN_Ca_old_) + (k_minus_ltrpn * LTRPN_Ca_old_)));	//14
    calc_i_CaL = (g_CaL * O_old_ * (V_old_ - E_CaL));	//22
    calc_i_pCa = ((i_pCa_max * pow(Cai_old_, 2.0e+00)) / (pow(Km_pCa, 2.0e+00) + pow(Cai_old_, 2.0e+00)));	//35
    calc_i_NaCa = (((((((k_NaCa * 1.0e+00) / (pow(K_mNa, 3.0e+00) + pow(Nao, 3.0e+00))) * 1.0e+00) / (K_mCa + Cao)) * 1.0e+00) / (1.0e+00 + (k_sat * exp((((eta - 1.0e+00) * V_old_ * F) / (R * T)))))) * ((exp(((eta * V_old_ * F) / (R * T))) * pow(Nai_old_, 3.0e+00) * Cao) - (exp((((eta - 1.0e+00) * V_old_ * F) / (R * T))) * pow(Nao, 3.0e+00) * Cai_old_)));	//36
    calc_E_CaN = (((R * T) / (2.0e+00 * F)) * log((Cao / Cai_old_)));	//38
    calc_E_Na = (((R * T) / F) * log((((9.0e-01 * Nao) + (1.0e-01 * Ko)) / ((9.0e-01 * Nai_old_) + (1.0e-01 * Ki_old_)))));	//41
    calc_E_K = (((R * T) / F) * log((Ko / Ki_old_)));	//68  
    calc_i_Kr = (g_Kr * O_K_old_ * (V_old_ - (((R * T) / F) * log((((9.80e-01 * Ko) + (2.0e-02 * Nao)) / ((9.80e-01 * Ki_old_) + (2.0e-02 * Nai_old_)))))));	//96
    calc_sigma = ((1.0e+00 / 7.0e+00) * (exp((Nao / 6.730e+04)) - 1.0e+00));	//110
    calc_O_ClCa = (2.0e-01 / (1.0e+00 + exp(((-(V_old_ - 4.670e+01)) / 7.80e+00))));	//112
    calc_i_Nab = (g_Nab * (V_old_ - calc_E_Na));	//65
    calc_i_Kto_s = (g_Kto_s * ato_s_old_ * ito_s_old_ * (V_old_ - calc_E_K));	//75
    calc_i_K1 = ((((2.9380e-01 * Ko) / (Ko + 2.10e+02)) * (V_old_ - calc_E_K)) / (1.0e+00 + exp((8.960e-02 * (V_old_ - calc_E_K)))));	//82
    calc_i_Ks = (g_Ks * pow(nKs_old_, 2.0e+00) * (V_old_ - calc_E_K));	//83
    calc_i_Kur = (g_Kur * aur_old_ * iur_old_ * (V_old_ - calc_E_K));	//87
    calc_i_Kss = (g_Kss * aKss_old_ * iKss_old_ * (V_old_ - calc_E_K));	//92
    calc_i_Cab = (g_Cab * (V_old_ - calc_E_CaN));	//37
    calc_i_Na = (g_Na * O_Na_old_ * (V_old_ - calc_E_Na));	//40
    calc_i_Kto_f = (g_Kto_f * pow(ato_f_old_, 3.0e+00) * ito_f_old_ * (V_old_ - calc_E_K));	//67
    calc_f_NaK = (1.0e+00 / (1.0e+00 + (1.2450e-01 * exp((((-1.0e-01) * V_old_ * F) / (R * T)))) + (3.650e-02 * calc_sigma * exp((((-V_old_) * F) / (R * T))))));	//109
    calc_i_ClCa = (((g_ClCa * calc_O_ClCa * Cai_old_) / (Cai_old_ + Km_Cl)) * (V_old_ - E_Cl));	//111
    calc_i_NaK = ((((i_NaK_max * calc_f_NaK * 1.0e+00) / (1.0e+00 + pow((Km_Nai / Nai_old_), 1.50e+00))) * Ko) / (Ko + Km_Ko));	//108 
}

void Bondarenko2004::calc_algs_hh(double* algs, double* pars, double* Y_old_, double time)
{
    calc_ass = (1.0e+00 / (1.0e+00 + exp(((-(V_old_ + 2.250e+01)) / 7.70e+00))));	//78
    calc_iss = (1.0e+00 / (1.0e+00 + exp(((V_old_ + 4.520e+01) / 5.70e+00))));	//79
    calc_tau_ta_s = ((4.930e-01 * exp(((-6.290e-02) * V_old_))) + 2.0580e+00);	//80
    calc_tau_ti_s = (2.70e+02 + (1.050e+03 / (1.0e+00 + exp(((V_old_ + 4.520e+01) / 5.70e+00)))));	//81
    calc_alpha_a = (1.80640e-01 * exp((3.5770e-02 * (V_old_ + 3.0e+01))));	//71
    calc_beta_a = (3.9560e-01 * exp(((-6.2370e-02) * (V_old_ + 3.0e+01))));	//72
    calc_alpha_i = ((1.520e-04 * exp(((-(V_old_ + 1.350e+01)) / 7.0e+00))) / ((6.70830e-02 * exp(((-(V_old_ + 3.350e+01)) / 7.0e+00))) + 1.0e+00));	//73
    calc_beta_i = ((9.50e-04 * exp(((V_old_ + 3.350e+01) / 7.0e+00))) / ((5.13350e-02 * exp(((V_old_ + 3.350e+01) / 7.0e+00))) + 1.0e+00));	//74
    calc_alpha_n = (abs(V_old_ + 2.650e+01) > EPSILON) ? 4.813330e-06 * (V_old_ + 2.650e+01) / (1.0e+00 - exp(-1.280e-01 * (V_old_ + 2.650e+01))) : 3.76040e-05;	//85
    calc_beta_n = (9.533330e-05 * exp(((-3.80e-02) * (V_old_ + 2.650e+01))));	//86
    calc_tau_aur = ((4.930e-01 * exp(((-6.290e-02) * V_old_))) + 2.0580e+00);	//90
    calc_tau_iur = (1.20e+03 - (1.70e+02 / (1.0e+00 + exp(((V_old_ + 4.520e+01) / 5.70e+00)))));	//91
    calc_tau_Kss = ((3.930e+01 * exp(((-8.620e-02) * V_old_))) + 1.3170e+01);	//95
}

void Bondarenko2004::calc_rhs_nl(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
    calc_algs_nl(algs, pars, Y_old_, t);

    V_f_ = ((-(calc_i_CaL + calc_i_pCa + calc_i_NaCa + calc_i_Cab + calc_i_Na + calc_i_Nab + calc_i_NaK + calc_i_Kto_f + calc_i_Kto_s + calc_i_K1 + calc_i_Ks + calc_i_Kur + calc_i_Kss + calc_i_Kr + calc_i_ClCa + calc_i_stim)));// 1
    Cai_f_ = ((calc_Bi * ((calc_J_leak + calc_J_xfer) - (calc_J_up + calc_J_trpn + ((((calc_i_Cab + calc_i_pCa) - (2.0e+00 * calc_i_NaCa)) * Acap * Cm) / (2.0e+00 * Vmyo * F))))));	// 2
    Cass_f_ = ((calc_Bss * (((calc_J_rel * VJSR) / Vss) - (((calc_J_xfer * Vmyo) / Vss) + ((calc_i_CaL * Acap * Cm) / (2.0e+00 * Vss * F))))));	// 3
    CaJSR_f_ = ((calc_BJSR * (calc_J_tr - calc_J_rel)));	// 4
    CaNSR_f_ = (((((calc_J_up - calc_J_leak) * Vmyo) / VNSR) - ((calc_J_tr * VJSR) / VNSR)));	// 5
    P_RyR_f_ = ((((-4.0e-02) * P_RyR_old_) - (((1.0e-01 * calc_i_CaL) / i_CaL_max) * exp(((-pow((V_old_ - 5.0e+00), 2.0e+00)) / 6.480e+02)))));	// 15
    LTRPN_Ca_f_ = (((k_plus_ltrpn * Cai_old_ * (LTRPN_tot - LTRPN_Ca_old_)) - (k_minus_ltrpn * LTRPN_Ca_old_)));	// 16
    HTRPN_Ca_f_ = (((k_plus_htrpn * Cai_old_ * (HTRPN_tot - HTRPN_Ca_old_)) - (k_minus_htrpn * HTRPN_Ca_old_)));	// 17
    Nai_f_ = ((((-(calc_i_Na + calc_i_Nab + (3.0e+00 * calc_i_NaK) + (3.0e+00 * calc_i_NaCa))) * Acap * Cm) / (Vmyo * F)));	// 39
    Ki_f_ = ((((-((calc_i_Kto_f + calc_i_Kto_s + calc_i_K1 + calc_i_Ks + calc_i_Kss + calc_i_Kur + calc_i_Kr) - (2.0e+00 * calc_i_NaK))) * Acap * Cm) / (Vmyo * F)));	// 66
    iKss_f_ = (0.0e+00);	// 94
}

void Bondarenko2004::calc_rhs_hh(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
    calc_algs_hh(algs, pars, Y_old_, t);

    ato_f_f_ = (((calc_alpha_a * (1.0e+00 - ato_f_old_)) - (calc_beta_a * ato_f_old_)));	// 69
    ito_f_f_ = (((calc_alpha_i * (1.0e+00 - ito_f_old_)) - (calc_beta_i * ito_f_old_)));	// 70
    ato_s_f_ = (((calc_ass - ato_s_old_) / calc_tau_ta_s));	// 76
    ito_s_f_ = (((calc_iss - ito_s_old_) / calc_tau_ti_s));	// 77
    nKs_f_ = (((calc_alpha_n * (1.0e+00 - nKs_old_)) - (calc_beta_n * nKs_old_)));	// 84
    aur_f_ = (((calc_ass - aur_old_) / calc_tau_aur));	// 88
    iur_f_ = (((calc_iss - iur_old_) / calc_tau_iur));	// 89
    aKss_f_ = (((calc_ass - aKss_old_) / calc_tau_Kss));	// 93
}

void Bondarenko2004::calc_hh_coeff(double* a, double* b, double* pars, double* algs, double* Y_old_, double t) 
{
    calc_algs_hh(algs, pars, Y_old_, t);

    ato_f_a_ = -calc_alpha_a - calc_beta_a;
    ito_f_a_ = -calc_alpha_i - calc_beta_i;	// 70
    ato_s_a_ = -1.0 / calc_tau_ta_s;	// 76
    ito_s_a_ = -1.0 / calc_tau_ti_s;	// 77
    nKs_a_ = -calc_alpha_n - calc_beta_n;	// 84
    aur_a_ = -1.0 / calc_tau_aur;	// 88
    iur_a_ = -1.0 / calc_tau_iur;	// 89
    aKss_a_ = -1.0 / calc_tau_Kss;	// 93

    ato_f_b_ = calc_alpha_a;	// 69
    ito_f_b_ = calc_alpha_i;	// 70
    ato_s_b_ = calc_ass / calc_tau_ta_s;	// 76
    ito_s_b_ = calc_iss / calc_tau_ti_s;	// 77
    nKs_b_ = calc_alpha_n;	// 84
    aur_b_ = calc_ass / calc_tau_aur;	// 88
    iur_b_ = calc_iss / calc_tau_iur;	// 89
    aKss_b_ = calc_ass / calc_tau_Kss;	// 93
}

void Bondarenko2004::prep_mk_transitions(double* algs, double* pars, double* Y_old_, double t, int mk_index) {
    if (mk_index == 0 || mk_index == -1) {
        calc_alpha_Na11 = (3.8020e+00 / ((1.0270e-01 * exp(((-(V_old_ + 2.50e+00)) / 1.70e+01))) + (2.0e-01 * exp(((-(V_old_ + 2.50e+00)) / 1.50e+02)))));	//51
        calc_alpha_Na12 = (3.8020e+00 / ((1.0270e-01 * exp(((-(V_old_ + 2.50e+00)) / 1.50e+01))) + (2.30e-01 * exp(((-(V_old_ + 2.50e+00)) / 1.50e+02)))));	//52
        calc_alpha_Na13 = (3.8020e+00 / ((1.0270e-01 * exp(((-(V_old_ + 2.50e+00)) / 1.20e+01))) + (2.50e-01 * exp(((-(V_old_ + 2.50e+00)) / 1.50e+02)))));	//53
        calc_beta_Na11 = (1.9170e-01 * exp(((-(V_old_ + 2.50e+00)) / 2.030e+01)));	//54
        calc_beta_Na12 = (2.0e-01 * exp(((-(V_old_ - 2.50e+00)) / 2.030e+01)));	//55
        calc_beta_Na13 = (2.20e-01 * exp(((-(V_old_ - 7.50e+00)) / 2.030e+01)));	//56
        calc_alpha_Na3 = (7.0e-07 * exp(((-(V_old_ + 7.0e+00)) / 7.70e+00)));	//57
        calc_beta_Na3 = (8.540000000000001e-03 + (2.0e-05 * V_old_));	//58
        calc_alpha_Na2 = (1.0e+00 / ((1.884950e-01 * exp(((-(V_old_ + 7.0e+00)) / 1.660e+01))) + 3.939560e-01));	//59
        calc_beta_Na2 = ((calc_alpha_Na13 * calc_alpha_Na2 * calc_alpha_Na3) / (calc_beta_Na13 * calc_beta_Na3));	//60
        calc_alpha_Na4 = (calc_alpha_Na2 / 1.0e+03);	//61
        calc_beta_Na4 = calc_alpha_Na3;	//62
        calc_alpha_Na5 = (calc_alpha_Na2 / 9.50e+04);	//63
        calc_beta_Na5 = (calc_alpha_Na3 / 5.0e+01);	//64
        C_Na3_old_ = (1.0e+00 - (O_Na_old_ + C_Na1_old_ + C_Na2_old_ + IF_Na_old_ + I1_Na_old_ + I2_Na_old_ + IC_Na2_old_ + IC_Na3_old_));
    }

    if (mk_index == 1 || mk_index == -1) {
        calc_alpha = ((4.0e-01 * exp(((V_old_ + 1.20e+01) / 1.0e+01)) * ((1.0e+00 + (7.0e-01 * exp(((-pow((V_old_ + 4.0e+01), 2.0e+00)) / 1.0e+01)))) - (7.50e-01 * exp(((-pow((V_old_ + 2.0e+01), 2.0e+00)) / 4.0e+02))))) / (1.0e+00 + (1.20e-01 * exp(((V_old_ + 1.20e+01) / 1.0e+01)))));	//31
        calc_beta = (5.0e-02 * exp(((-(V_old_ + 1.20e+01)) / 1.30e+01)));	//32
        calc_gamma = ((Kpc_max * Cass_old_) / (Kpc_half + Cass_old_));	//33
        calc_Kpcf = (1.30e+01 * (1.0e+00 - exp(((-pow((V_old_ + 1.450e+01), 2.0e+00)) / 1.0e+02))));	//34
        C1_old_ = (1.0e+00 - (O_old_ + C2_old_ + C3_old_ + C4_old_ + I1_old_ + I2_old_ + I3_old_));	//24
    }

    if (mk_index == 2 || mk_index == -1) {
        calc_alpha_a0 = (2.23480e-02 * exp((1.1760e-02 * V_old_)));	//102
        calc_beta_a0 = (4.70020e-02 * exp(((-6.310e-02) * V_old_)));	//103
        calc_alpha_a1 = (1.37330e-02 * exp((3.81980e-02 * V_old_)));	//104
        calc_beta_a1 = (6.889999999999999e-05 * exp(((-4.1780e-02) * V_old_)));	//105
        calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current = (9.08210e-02 * exp((2.33910e-02 * (V_old_ + 5.0e+00))));	//106
        calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current = (6.4970e-03 * exp(((-3.2680e-02) * (V_old_ + 5.0e+00))));	//107
        C_K0_old_ = (1.0e+00 - (C_K1_old_ + C_K2_old_ + O_K_old_ + I_K_old_));
    }

    if (mk_index == 3 || mk_index == -1) {
        P_C1_old_ = (1.0e+00 - (P_C2_old_ + P_O1_old_ + P_O2_old_));	//19
    }
}

void Bondarenko2004::calc_rhs_mk(double* rhs, double* pars, double* algs, double* Y_old_, double t, int mk_index)
{
    prep_mk_transitions(algs, pars, Y_old_, t, mk_index);

    if (mk_index == 0 || mk_index == -1) {
        C_Na2_f_ = ((((calc_alpha_Na11 * C_Na3_old_) + (calc_beta_Na12 * C_Na1_old_) + (calc_alpha_Na3 * IC_Na2_old_)) - ((calc_beta_Na11 * C_Na2_old_) + (calc_alpha_Na12 * C_Na2_old_) + (calc_beta_Na3 * C_Na2_old_))));	// 43
        C_Na1_f_ = ((((calc_alpha_Na12 * C_Na2_old_) + (calc_beta_Na13 * O_Na_old_) + (calc_alpha_Na3 * IF_Na_old_)) - ((calc_beta_Na12 * C_Na1_old_) + (calc_alpha_Na13 * C_Na1_old_) + (calc_beta_Na3 * C_Na1_old_))));	// 44
        O_Na_f_ = ((((calc_alpha_Na13 * C_Na1_old_) + (calc_beta_Na2 * IF_Na_old_)) - ((calc_beta_Na13 * O_Na_old_) + (calc_alpha_Na2 * O_Na_old_))));	// 45
        IF_Na_f_ = ((((calc_alpha_Na2 * O_Na_old_) + (calc_beta_Na3 * C_Na1_old_) + (calc_beta_Na4 * I1_Na_old_) + (calc_alpha_Na12 * IC_Na2_old_)) - ((calc_beta_Na2 * IF_Na_old_) + (calc_alpha_Na3 * IF_Na_old_) + (calc_alpha_Na4 * IF_Na_old_) + (calc_beta_Na12 * IF_Na_old_))));	// 46
        I1_Na_f_ = ((((calc_alpha_Na4 * IF_Na_old_) + (calc_beta_Na5 * I2_Na_old_)) - ((calc_beta_Na4 * I1_Na_old_) + (calc_alpha_Na5 * I1_Na_old_))));	// 47
        I2_Na_f_ = (((calc_alpha_Na5 * I1_Na_old_) - (calc_beta_Na5 * I2_Na_old_)));	// 48
        IC_Na2_f_ = ((((calc_alpha_Na11 * IC_Na3_old_) + (calc_beta_Na12 * IF_Na_old_) + (calc_beta_Na3 * C_Na2_old_)) - ((calc_beta_Na11 * IC_Na2_old_) + (calc_alpha_Na12 * IC_Na2_old_) + (calc_alpha_Na3 * IC_Na2_old_))));	// 49
        IC_Na3_f_ = ((((calc_beta_Na11 * IC_Na2_old_) + (calc_beta_Na3 * C_Na3_old_)) - ((calc_alpha_Na11 * IC_Na3_old_) + (calc_alpha_Na3 * IC_Na3_old_))));	// 50
        C_Na3_f_ = (1.0e+00 - (O_Na_f_ + C_Na1_f_ + C_Na2_f_ + IF_Na_f_ + I1_Na_f_ + I2_Na_f_ + IC_Na2_f_ + IC_Na3_f_));
    }

    if (mk_index == 1 || mk_index == -1) {
        O_f_ = ((((calc_alpha * C4_old_) + (Kpcb * I1_old_) + (1.0e-03 * ((calc_alpha * I2_old_) - (calc_Kpcf * O_old_)))) - ((4.0e+00 * calc_beta * O_old_) + (calc_gamma * O_old_))));	// 23
        C2_f_ = ((((4.0e+00 * calc_alpha * C1_old_) + (2.0e+00 * calc_beta * C3_old_)) - ((calc_beta * C2_old_) + (3.0e+00 * calc_alpha * C2_old_))));	// 25
        C3_f_ = ((((3.0e+00 * calc_alpha * C2_old_) + (3.0e+00 * calc_beta * C4_old_)) - ((2.0e+00 * calc_beta * C3_old_) + (2.0e+00 * calc_alpha * C3_old_))));	// 26
        C4_f_ = ((((2.0e+00 * calc_alpha * C3_old_) + (4.0e+00 * calc_beta * O_old_) + (1.0e-02 * ((4.0e+00 * Kpcb * calc_beta * I1_old_) - (calc_alpha * calc_gamma * C4_old_))) + (2.0e-03 * ((4.0e+00 * calc_beta * I2_old_) - (calc_Kpcf * C4_old_))) + (4.0e+00 * calc_beta * Kpcb * I3_old_)) - ((3.0e+00 * calc_beta * C4_old_) + (calc_alpha * C4_old_) + (1.0e+00 * calc_gamma * calc_Kpcf * C4_old_))));	// 27
        I1_f_ = ((((calc_gamma * O_old_) + (1.0e-03 * ((calc_alpha * I3_old_) - (calc_Kpcf * I1_old_))) + (1.0e-02 * ((calc_alpha * calc_gamma * C4_old_) - (4.0e+00 * calc_beta * Kpcb * I1_old_)))) - (Kpcb * I1_old_)));	// 28
        I2_f_ = ((((1.0e-03 * ((calc_Kpcf * O_old_) - (calc_alpha * I2_old_))) + (Kpcb * I3_old_) + (2.0e-03 * ((calc_Kpcf * C4_old_) - (4.0e+00 * calc_beta * I2_old_)))) - (calc_gamma * I2_old_)));	// 29
        I3_f_ = ((((1.0e-03 * ((calc_Kpcf * I1_old_) - (calc_alpha * I3_old_))) + (calc_gamma * I2_old_) + (1.0e+00 * calc_gamma * calc_Kpcf * C4_old_)) - ((4.0e+00 * calc_beta * Kpcb * I3_old_) + (Kpcb * I3_old_))));	// 30
        C1_f_ = (1.0e+00 - (O_f_ + C2_f_ + C3_f_ + C4_f_ + I1_f_ + I2_f_ + I3_f_));	//24
    }

    if (mk_index == 2 || mk_index == -1) {
        C_K2_f_ = ((((kf * C_K1_old_) + (calc_beta_a1 * O_K_old_)) - ((kb * C_K2_old_) + (calc_alpha_a1 * C_K2_old_))));	// 98
        C_K1_f_ = ((((calc_alpha_a0 * C_K0_old_) + (kb * C_K2_old_)) - ((calc_beta_a0 * C_K1_old_) + (kf * C_K1_old_))));	// 99
        O_K_f_ = ((((calc_alpha_a1 * C_K2_old_) + (calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current * I_K_old_)) - ((calc_beta_a1 * O_K_old_) + (calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current * O_K_old_))));	// 100
        I_K_f_ = (((calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current * O_K_old_) - (calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current * I_K_old_)));
        C_K0_f_ = (1.0e+00 - (C_K1_f_ + C_K2_f_ + O_K_f_ + I_K_f_));
    }

    if (mk_index == 3 || mk_index == -1) {
        P_O1_f_ = ((((k_plus_a * pow(Cass_old_, n) * P_C1_old_) + (k_minus_b * P_O2_old_) + (k_minus_c * P_C2_old_)) - ((k_minus_a * P_O1_old_) + (k_plus_b * pow(Cass_old_, m) * P_O1_old_) + (k_plus_c * P_O1_old_))));	// 18
        P_O2_f_ = (((k_plus_b * pow(Cass_old_, m) * P_O1_old_) - (k_minus_b * P_O2_old_)));	// 20
        P_C2_f_ = (((k_plus_c * P_O1_old_) - (k_minus_c * P_C2_old_)));	// 21
        P_C1_f_ = (1.0e+00 - (P_C2_f_ + P_O1_f_ + P_O2_f_));	//19
    }
}

void Bondarenko2004::calc_mk_transitions(double** Tr, int mk_index, double* pars, double* algs, double* Y_old_, double t) 
{
    switch (mk_index) {
        case 0:
            //C_Na3
            Tr[0][0] = 0;
            Tr[1][0] = calc_alpha_Na11;
            Tr[8][0] = calc_beta_Na3;

            //C_Na2
            Tr[0][1] = calc_beta_Na11;
            Tr[1][1] = 0;
            Tr[2][1] = calc_alpha_Na12;
            Tr[7][1] = calc_beta_Na3;

            //C_Na1
            Tr[1][2] = calc_beta_Na12;
            Tr[2][2] = 0;
            Tr[3][2] = calc_alpha_Na13;
            Tr[4][2] = calc_beta_Na3;

            //O_Na
            Tr[2][3] = calc_beta_Na13;
            Tr[3][3] = 0;
            Tr[4][3] = calc_alpha_Na2;

            //IF_Na
            Tr[2][4] = calc_alpha_Na3;
            Tr[3][4] = calc_beta_Na2;
            Tr[4][4] = 0;
            Tr[5][4] = calc_alpha_Na4;
            Tr[7][4] = calc_beta_Na12;

            //I1_Na
            Tr[4][5] = calc_beta_Na4;
            Tr[5][5] = 0;
            Tr[6][5] = calc_alpha_Na5;

            //I2_Na
            Tr[5][6] = calc_beta_Na5;
            Tr[6][6] = 0;

            //IC_Na2
            Tr[1][7] = calc_alpha_Na3;
            Tr[4][7] = calc_alpha_Na12;
            Tr[7][7] = 0;
            Tr[8][7] = calc_beta_Na11;

            //IC_Na3
            Tr[0][8] = calc_alpha_Na3;
            Tr[7][8] = calc_alpha_Na11;
            Tr[8][8] = 0;
            break;
        case 1:
            Tr[1][0] = 4 * calc_alpha;

            Tr[0][1] = calc_beta;
            Tr[2][1] = 3 * calc_alpha;
            //C3
            Tr[1][2] = 3 * calc_beta;
            Tr[3][2] = 2 * calc_alpha;
            //C4
            Tr[2][3] = 3 * calc_beta;
            Tr[4][3] = 0.01 * (calc_alpha * calc_gamma);
            Tr[5][3] = 0.002 * calc_Kpcf;
            Tr[6][3] = calc_gamma * calc_Kpcf;
            Tr[7][3] = calc_alpha;
            //I1
            Tr[3][4] = 0.01 * (4 * Kpcb * calc_beta);
            Tr[6][4] = 0.001 * calc_Kpcf;
            Tr[7][4] = Kpcb;
            //I2
            Tr[3][5] = 0.002 * (4 * calc_beta);
            Tr[6][5] = calc_gamma;
            Tr[7][5] = 0.001 * calc_alpha;
            //I3
            Tr[3][6] = 4 * calc_beta * Kpcb;
            Tr[4][6] = 0.001 * calc_alpha;
            Tr[5][6] = Kpcb;
            //O
            Tr[3][7] = 4 * calc_beta;
            Tr[4][7] = calc_gamma;
            Tr[5][7] = 0.001 * calc_Kpcf;
            break;
        case 2:
            Tr[1][0] = (calc_alpha_a0);

            Tr[0][1] = calc_beta_a0;
            Tr[2][1] = kf;

            Tr[1][2] = kb;
            Tr[3][2] = calc_alpha_a1;

            Tr[2][3] = calc_beta_a1;
            Tr[4][3] = calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current;

            Tr[3][4] = calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current;
            break;
        case 3:
            //P_C1d
            Tr[1][0] = k_plus_a * pow(Cass_old_, n);
            //P_O1d
            Tr[0][1] = k_minus_a;
            Tr[2][1] = k_plus_b * pow(Cass_old_, m);
            Tr[3][1] = k_plus_c;
            //P_O2d
            Tr[1][2] = k_minus_b;
            //P_C2d
            Tr[1][3] = k_minus_c;
            break;
    }
}

bool Bondarenko2004::has_single_rhs_formula(int i) {
    return i == 0 || i == 2;
}
	
double Bondarenko2004::calc_single_rhs_formula(int i, double* pars, double* Y_old_, double t) {
    if (i == 0) {
        double c_i_stim = calc_stimulus(pars, t);
        double c_i_CaL = (g_CaL*O_old_*(V_old_-E_CaL));
        double c_i_pCa = ((i_pCa_max*pow(Cai_old_,2.000000000000000e+00))/(pow(Km_pCa,2.000000000000000e+00)+pow(Cai_old_,2.000000000000000e+00)));
        double c_i_NaCa = (((((((k_NaCa*1.000000000000000e+00)/(pow(K_mNa,3.000000000000000e+00)+pow(Nao,3.000000000000000e+00)))*1.000000000000000e+00)/(K_mCa+Cao))*1.000000000000000e+00)/(1.000000000000000e+00+(k_sat*exp((((eta-1.000000000000000e+00)*V_old_*F)/(R*T))))))*((exp(((eta*V_old_*F)/(R*T)))*pow(Nai_old_,3.000000000000000e+00)*Cao)-(exp((((eta-1.000000000000000e+00)*V_old_*F)/(R*T)))*pow(Nao,3.000000000000000e+00)*Cai_old_)));
        double c_E_CaN = (((R*T)/(2.000000000000000e+00*F))*log((Cao/Cai_old_)));
        double c_E_Na = (((R*T)/F)*log((((9.000000000000000e-01*Nao)+(1.000000000000000e-01*Ko))/((9.000000000000000e-01*Nai_old_)+(1.000000000000000e-01*Ki_old_)))));
        double c_E_K = (((R*T)/F)*log((Ko/Ki_old_)));
        double c_i_Kr = (g_Kr*O_K_old_*(V_old_-(((R*T)/F)*log((((9.800000000000000e-01*Ko)+(2.000000000000000e-02*Nao))/((9.800000000000000e-01*Ki_old_)+(2.000000000000000e-02*Nai_old_)))))));
        double c_sigma = ((1.000000000000000e+00/7.000000000000000e+00)*(exp((Nao/6.730000000000000e+04))-1.000000000000000e+00));
        double c_O_ClCa = (2.000000000000000e-01/(1.000000000000000e+00+exp(((-(V_old_-4.670000000000000e+01))/7.800000000000000e+00))));
        double c_i_Nab = (g_Nab*(V_old_-c_E_Na));
        double c_i_Kto_s = (g_Kto_s*ato_s_old_*ito_s_old_*(V_old_-c_E_K));
        double c_i_K1 = ((((2.938000000000000e-01*Ko)/(Ko+2.100000000000000e+02))*(V_old_-c_E_K))/(1.000000000000000e+00+exp((8.960000000000000e-02*(V_old_-c_E_K)))));
        double c_i_Ks = (g_Ks*pow(nKs_old_,2.000000000000000e+00)*(V_old_-c_E_K));
        double c_i_Kur = (g_Kur*aur_old_*iur_old_*(V_old_-c_E_K));
        double c_i_Kss = (g_Kss*aKss_old_*iKss_old_*(V_old_-c_E_K));
        double c_i_Cab = (g_Cab*(V_old_-c_E_CaN));
        double c_i_Na = (g_Na*O_Na_old_*(V_old_-c_E_Na));
        double c_i_Kto_f = (g_Kto_f*pow(ato_f_old_,3.000000000000000e+00)*ito_f_old_*(V_old_-c_E_K));
        double c_f_NaK = (1.000000000000000e+00/(1.000000000000000e+00+(1.245000000000000e-01*exp((((-1.000000000000000e-01)*V_old_*F)/(R*T))))+(3.650000000000000e-02*c_sigma*exp((((-V_old_)*F)/(R*T))))));
        double c_i_ClCa = (((g_ClCa*c_O_ClCa*Cai_old_)/(Cai_old_+Km_Cl))*(V_old_-E_Cl));
        double c_i_NaK = ((((i_NaK_max*c_f_NaK*1.000000000000000e+00)/(1.000000000000000e+00+pow((Km_Nai/Nai_old_),1.500000000000000e+00)))*Ko)/(Ko+Km_Ko));
        return ((-(c_i_CaL+c_i_pCa+c_i_NaCa+c_i_Cab+c_i_Na+c_i_Nab+c_i_NaK+c_i_Kto_f+c_i_Kto_s+c_i_K1+c_i_Ks+c_i_Kur+c_i_Kss+c_i_Kr+c_i_ClCa+c_i_stim)));
    } 
    /*
    else if (i == 1) {
        double c_Bi = pow((1.000000000000000e+00+((CMDN_tot*Km_CMDN)/pow((Km_CMDN+Cai_old_),2.000000000000000e+00))),(-1.000000000000000e+00));
        double c_J_xfer = ((Cass_old_-Cai_old_)/tau_xfer);
        double c_J_leak = (v2*(CaNSR_old_-Cai_old_));
        double c_J_up = ((v3*pow(Cai_old_,2.000000000000000e+00))/(pow(Km_up,2.000000000000000e+00)+pow(Cai_old_,2.000000000000000e+00)));
        double c_J_trpn = (((k_plus_htrpn*Cai_old_*(HTRPN_tot-HTRPN_Ca_old_))+(k_plus_ltrpn*Cai_old_*(LTRPN_tot-LTRPN_Ca_old_)))-((k_minus_htrpn*HTRPN_Ca_old_)+(k_minus_ltrpn*LTRPN_Ca_old_)));
        double c_i_pCa = ((i_pCa_max*pow(Cai_old_,2.000000000000000e+00))/(pow(Km_pCa,2.000000000000000e+00)+pow(Cai_old_,2.000000000000000e+00)));
        double c_i_NaCa = (((((((k_NaCa*1.000000000000000e+00)/(pow(K_mNa,3.000000000000000e+00)+pow(Nao,3.000000000000000e+00)))*1.000000000000000e+00)/(K_mCa+Cao))*1.000000000000000e+00)/(1.000000000000000e+00+(k_sat*exp((((eta-1.000000000000000e+00)*V_old_*F)/(R*T))))))*((exp(((eta*V_old_*F)/(R*T)))*pow(Nai_old_,3.000000000000000e+00)*Cao)-(exp((((eta-1.000000000000000e+00)*V_old_*F)/(R*T)))*pow(Nao,3.000000000000000e+00)*Cai_old_)));
        double c_E_CaN = (((R*T)/(2.000000000000000e+00*F))*log((Cao/Cai_old_)));
        double c_i_Cab = (g_Cab*(V_old_-c_E_CaN));
        return ((c_Bi*((c_J_leak+c_J_xfer)-(c_J_up+c_J_trpn+((((c_i_Cab+c_i_pCa)-(2.000000000000000e+00*c_i_NaCa))*Acap*Cm)/(2.000000000000000e+00*Vmyo*F))))));
    } /* */
    else if (i == 2) {
        double c_Bss = pow((1.0+((CMDN_tot*Km_CMDN)/pow((Km_CMDN+Cass_old_),2.0))),(-1.0));
        double c_J_rel = (v1*(P_O1_old_+P_O2_old_)*(CaJSR_old_-Cass_old_)*P_RyR_old_);
        double c_J_xfer = ((Cass_old_-Cai_old_)/tau_xfer);
        double c_i_CaL = (g_CaL*O_old_*(V_old_-E_CaL));
        return ((c_Bss*(((c_J_rel*VJSR)/Vss)-(((c_J_xfer*Vmyo)/Vss)+((c_i_CaL*Acap*Cm)/(2.0*Vss*F))))));
    }

    return 0;
}

bool Bondarenko2004::is_mk_active(int i) {
    return i == 3;
}