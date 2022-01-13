#pragma once

#define stim_state pars[0]
#define stim_amplitude pars[1] 	 // picoA_per_picoF
#define stim_start pars[2] 	 // millisecond
#define stim_end pars[3] 	 // millisecond
#define stim_period pars[4] 	 // millisecond
#define stim_duration pars[5] 	 // millisecond
#define Acap pars[6] 	 // cm2
#define Cm pars[7] 	 // microF_per_cm2
#define Vmyo pars[8] 	 // microlitre
#define F pars[9] 	 // coulomb_per_millimole
#define VJSR pars[10] 	 // microlitre
#define Vss pars[11] 	 // microlitre
#define VNSR pars[12] 	 // microlitre
#define CMDN_tot pars[13] 	 // micromolar
#define Km_CMDN pars[14] 	 // micromolar
#define CSQN_tot pars[15] 	 // micromolar
#define Km_CSQN pars[16] 	 // micromolar
#define v1 pars[17] 	 // per_millisecond
#define tau_tr pars[18] 	 // millisecond
#define tau_xfer pars[19] 	 // millisecond
#define v2 pars[20] 	 // per_millisecond
#define v3 pars[21] 	 // micromolar_per_millisecond
#define Km_up pars[22] 	 // micromolar
#define k_plus_htrpn pars[23] 	 // per_micromolar_millisecond
#define HTRPN_tot pars[24] 	 // micromolar
#define k_plus_ltrpn pars[25] 	 // per_micromolar_millisecond
#define LTRPN_tot pars[26] 	 // micromolar
#define k_minus_htrpn pars[27] 	 // per_millisecond
#define k_minus_ltrpn pars[28] 	 // per_millisecond
#define i_CaL_max pars[29] 	 // picoA_per_picoF
#define k_plus_a pars[30] 	 // micromolar4_per_millisecond
#define n pars[31] 	 // dimensionless
#define k_minus_b pars[32] 	 // per_millisecond
#define k_minus_c pars[33] 	 // per_millisecond
#define k_minus_a pars[34] 	 // per_millisecond
#define k_plus_b pars[35] 	 // micromolar3_per_millisecond
#define m pars[36] 	 // dimensionless
#define k_plus_c pars[37] 	 // per_millisecond
#define g_CaL pars[38] 	 // milliS_per_microF
#define E_CaL pars[39] 	 // millivolt
#define Kpcb pars[40] 	 // per_millisecond
#define Kpc_max pars[41] 	 // per_millisecond
#define Kpc_half pars[42] 	 // micromolar
#define i_pCa_max pars[43] 	 // picoA_per_picoF
#define Km_pCa pars[44] 	 // micromolar
#define k_NaCa pars[45] 	 // picoA_per_picoF
#define K_mNa pars[46] 	 // micromolar
#define Nao pars[47] 	 // micromolar
#define K_mCa pars[48] 	 // micromolar
#define Cao pars[49] 	 // micromolar
#define k_sat pars[50] 	 // dimensionless
#define eta pars[51] 	 // dimensionless
#define R pars[52] 	 // joule_per_mole_kelvin
#define T pars[53] 	 // kelvin
#define g_Cab pars[54] 	 // milliS_per_microF
#define g_Na pars[55] 	 // milliS_per_microF
#define Ko pars[56] 	 // micromolar
#define g_Nab pars[57] 	 // milliS_per_microF
#define g_Kto_f pars[58] 	 // milliS_per_microF
#define g_Kto_s pars[59] 	 // milliS_per_microF
#define g_Ks pars[60] 	 // milliS_per_microF
#define g_Kur pars[61] 	 // milliS_per_microF
#define g_Kss pars[62] 	 // milliS_per_microF
#define g_Kr pars[63] 	 // milliS_per_microF
#define kf pars[64] 	 // per_millisecond
#define kb pars[65] 	 // per_millisecond
#define i_NaK_max pars[66] 	 // picoA_per_picoF
#define Km_Nai pars[67] 	 // micromolar
#define Km_Ko pars[68] 	 // micromolar
#define g_ClCa pars[69] 	 // milliS_per_microF
#define Km_Cl pars[70] 	 // micromolar
#define E_Cl pars[71] 	 // millivolt

#define calc_i_stim algs[0] 	 // picoA_per_picoF
#define calc_Bi algs[1] 	 // dimensionless
#define calc_Bss algs[2] 	 // dimensionless
#define calc_BJSR algs[3] 	 // dimensionless
#define calc_J_rel algs[4] 	 // micromolar_per_millisecond
#define calc_J_tr algs[5] 	 // micromolar_per_millisecond
#define calc_J_xfer algs[6] 	 // micromolar_per_millisecond
#define calc_J_leak algs[7] 	 // micromolar_per_millisecond
#define calc_J_up algs[8] 	 // micromolar_per_millisecond
#define calc_J_trpn algs[9] 	 // micromolar_per_millisecond
#define calc_i_CaL algs[10] 	 // picoA_per_picoF
#define calc_alpha algs[11] 	 // per_millisecond
#define calc_beta algs[12] 	 // per_millisecond
#define calc_gamma algs[13] 	 // per_millisecond
#define calc_Kpcf algs[14] 	 // per_millisecond
#define calc_i_pCa algs[15] 	 // picoA_per_picoF
#define calc_i_NaCa algs[16] 	 // picoA_per_picoF
#define calc_i_Cab algs[17] 	 // picoA_per_picoF
#define calc_E_CaN algs[18] 	 // millivolt
#define calc_i_Na algs[19] 	 // picoA_per_picoF
#define calc_E_Na algs[20] 	 // millivolt
#define calc_alpha_Na11 algs[21] 	 // per_millisecond
#define calc_alpha_Na12 algs[22] 	 // per_millisecond
#define calc_alpha_Na13 algs[23] 	 // per_millisecond
#define calc_beta_Na11 algs[24] 	 // per_millisecond
#define calc_beta_Na12 algs[25] 	 // per_millisecond
#define calc_beta_Na13 algs[26] 	 // per_millisecond
#define calc_alpha_Na3 algs[27] 	 // per_millisecond
#define calc_beta_Na3 algs[28] 	 // per_millisecond
#define calc_alpha_Na2 algs[29] 	 // per_millisecond
#define calc_beta_Na2 algs[30] 	 // per_millisecond
#define calc_alpha_Na4 algs[31] 	 // per_millisecond
#define calc_beta_Na4 algs[32] 	 // per_millisecond
#define calc_alpha_Na5 algs[33] 	 // per_millisecond
#define calc_beta_Na5 algs[34] 	 // per_millisecond
#define calc_i_Nab algs[35] 	 // picoA_per_picoF
#define calc_i_Kto_f algs[36] 	 // picoA_per_picoF
#define calc_E_K algs[37] 	 // millivolt
#define calc_alpha_a algs[38] 	 // per_millisecond
#define calc_beta_a algs[39] 	 // per_millisecond
#define calc_alpha_i algs[40] 	 // per_millisecond
#define calc_beta_i algs[41] 	 // per_millisecond
#define calc_i_Kto_s algs[42] 	 // picoA_per_picoF
#define calc_ass algs[43] 	 // dimensionless
#define calc_iss algs[44] 	 // dimensionless
#define calc_tau_ta_s algs[45] 	 // millisecond
#define calc_tau_ti_s algs[46] 	 // millisecond
#define calc_i_K1 algs[47] 	 // picoA_per_picoF
#define calc_i_Ks algs[48] 	 // picoA_per_picoF
#define calc_alpha_n algs[49] 	 // per_millisecond
#define calc_beta_n algs[50] 	 // per_millisecond
#define calc_i_Kur algs[51] 	 // picoA_per_picoF
#define calc_tau_aur algs[52] 	 // millisecond
#define calc_tau_iur algs[53] 	 // millisecond
#define calc_i_Kss algs[54] 	 // picoA_per_picoF
#define calc_tau_Kss algs[55] 	 // millisecond
#define calc_i_Kr algs[56] 	 // picoA_per_picoF
#define calc_alpha_a0 algs[57] 	 // per_millisecond
#define calc_beta_a0 algs[58] 	 // per_millisecond
#define calc_alpha_a1 algs[59] 	 // per_millisecond
#define calc_beta_a1 algs[60] 	 // per_millisecond
#define calc_i_NaK algs[61] 	 // picoA_per_picoF
#define calc_f_NaK algs[62] 	 // dimensionless
#define calc_sigma algs[63] 	 // dimensionless
#define calc_i_ClCa algs[64] 	 // picoA_per_picoF
#define calc_O_ClCa algs[65] 	 // dimensionless
#define calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current algs[66] 	 // (null)
#define calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current algs[67] 	 // (null)

#define V_old_ Y_old_[0]
#define Cai_old_ Y_old_[1]
#define Cass_old_ Y_old_[2]
#define CaJSR_old_ Y_old_[3]
#define CaNSR_old_ Y_old_[4]
#define P_RyR_old_ Y_old_[5]
#define LTRPN_Ca_old_ Y_old_[6]
#define HTRPN_Ca_old_ Y_old_[7]
#define Nai_old_ Y_old_[8]
#define Ki_old_ Y_old_[9]
#define iKss_old_ Y_old_[10]
#define ato_f_old_ Y_old_[11]
#define ito_f_old_ Y_old_[12]
#define ato_s_old_ Y_old_[13]
#define ito_s_old_ Y_old_[14]
#define nKs_old_ Y_old_[15]
#define aur_old_ Y_old_[16]
#define iur_old_ Y_old_[17]
#define aKss_old_ Y_old_[18]
#define C_Na3_old_ Y_old_[19]
#define C_Na2_old_ Y_old_[20]
#define C_Na1_old_ Y_old_[21]
#define O_Na_old_ Y_old_[22]
#define IF_Na_old_ Y_old_[23]
#define I1_Na_old_ Y_old_[24]
#define I2_Na_old_ Y_old_[25]
#define IC_Na2_old_ Y_old_[26]
#define IC_Na3_old_ Y_old_[27]
#define C1_old_ Y_old_[28]
#define C2_old_ Y_old_[29]
#define C3_old_ Y_old_[30]
#define C4_old_ Y_old_[31]
#define I1_old_ Y_old_[32]
#define I2_old_ Y_old_[33]
#define I3_old_ Y_old_[34]
#define O_old_ Y_old_[35]
#define C_K0_old_ Y_old_[36]
#define C_K2_old_ Y_old_[37]
#define C_K1_old_ Y_old_[38]
#define O_K_old_ Y_old_[39]
#define I_K_old_ Y_old_[40]
#define P_C1_old_ Y_old_[41]
#define P_O1_old_ Y_old_[42]
#define P_O2_old_ Y_old_[43]
#define P_C2_old_ Y_old_[44]

#define V_f_ rhs[0]
#define Cai_f_ rhs[1]
#define Cass_f_ rhs[2]
#define CaJSR_f_ rhs[3]
#define CaNSR_f_ rhs[4]
#define P_RyR_f_ rhs[5]
#define LTRPN_Ca_f_ rhs[6]
#define HTRPN_Ca_f_ rhs[7]
#define Nai_f_ rhs[8]
#define Ki_f_ rhs[9]
#define iKss_f_ rhs[10]
#define ato_f_f_ rhs[11]
#define ito_f_f_ rhs[12]
#define ato_s_f_ rhs[13]
#define ito_s_f_ rhs[14]
#define nKs_f_ rhs[15]
#define aur_f_ rhs[16]
#define iur_f_ rhs[17]
#define aKss_f_ rhs[18]
#define C_Na3_f_ rhs[19]
#define C_Na2_f_ rhs[20]
#define C_Na1_f_ rhs[21]
#define O_Na_f_ rhs[22]
#define IF_Na_f_ rhs[23]
#define I1_Na_f_ rhs[24]
#define I2_Na_f_ rhs[25]
#define IC_Na2_f_ rhs[26]
#define IC_Na3_f_ rhs[27]
#define C1_f_ rhs[28]
#define C2_f_ rhs[29]
#define C3_f_ rhs[30]
#define C4_f_ rhs[31]
#define I1_f_ rhs[32]
#define I2_f_ rhs[33]
#define I3_f_ rhs[34]
#define O_f_ rhs[35]
#define C_K0_f_ rhs[36]
#define C_K2_f_ rhs[37]
#define C_K1_f_ rhs[38]
#define O_K_f_ rhs[39]
#define I_K_f_ rhs[40]
#define P_C1_f_ rhs[41]
#define P_O1_f_ rhs[42]
#define P_O2_f_ rhs[43]
#define P_C2_f_ rhs[44]

#define ato_f_a_ a[0]
#define ito_f_a_ a[1]
#define ato_s_a_ a[2]
#define ito_s_a_ a[3]
#define nKs_a_ a[4]
#define aur_a_ a[5]
#define iur_a_ a[6]
#define aKss_a_ a[7]

#define ato_f_b_ b[0]
#define ito_f_b_ b[1]
#define ato_s_b_ b[2]
#define ito_s_b_ b[3]
#define nKs_b_ b[4]
#define aur_b_ b[5]
#define iur_b_ b[6]
#define aKss_b_ b[7]