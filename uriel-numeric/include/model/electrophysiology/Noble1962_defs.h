#pragma once

#define stim_state pars[0] // This must always be the first param of a myocyte model
#define stim_amplitude pars[1] // And this must always be the second
#define stim_period pars[2]
#define stim_start pars[3]
#define stim_duration pars[4]
#define Cm pars[5]
#define g_Na_max pars[6]
#define E_Na pars[7]
#define g_L pars[8]
#define E_L pars[9]

#define V_old_ Y_old_[0]
#define m_old_ Y_old_[1]
#define h_old_ Y_old_[2]
#define n_old_ Y_old_[3]

#define V_f_ rhs[0]
#define m_f_ rhs[1]
#define h_f_ rhs[2]
#define n_f_ rhs[3]

#define m_a_ a[0] 
#define h_a_ a[1] 
#define n_a_ a[2]

#define m_b_ b[0] 
#define h_b_ b[1] 
#define n_b_ b[2]

#define calc_i_Stim algs[0] 	 
#define calc_g_Na algs[1] 	 
#define calc_alpha_m algs[2]
#define calc_beta_m algs[3]
#define calc_m_inf algs[4]
#define calc_tau_m algs[5]
#define calc_alpha_h algs[6]
#define calc_beta_h algs[7]
#define calc_h_inf algs[8]
#define calc_tau_h algs[9]
#define calc_alpha_n algs[10]
#define calc_beta_n algs[11]
#define calc_n_inf algs[12]
#define calc_tau_n algs[13]
#define calc_g_K1 algs[14]
#define calc_g_K2 algs[15]
#define calc_i_Na algs[16]
#define calc_i_K algs[17]
#define calc_i_Leak algs[18]