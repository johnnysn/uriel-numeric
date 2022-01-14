#include "pch.h"
#include "model/electrophysiology/ToRORd.h"
#include "model/electrophysiology/ToRORd_defs.h"

ToRORd::ToRORd() : CellModel(11, 8, 26, 1, 68, 72)
{
    nStates_MKM[0] = 9;

    nStates_MKM_max = 9;
}

void ToRORd::set_default_parameters(double* pars)
{
    celltype = 0;
    nao = 140.0;
    cao = 1.8;
    ko = 5.0;
    clo = 150.0;
    R = 8314;
    T = 310;
    F = 96485;
    zna = 1;
    zca = 2;
    zk = 1;
    zcl = -1;
    L = 0.01;
    rad = 0.0011;
    i_Stim_Start = 0;
    i_Stim_End = 100000000000000000;
    i_Stim_Amplitude = -53;
    i_Stim_Period = 1000;
    i_Stim_PulseDuration = 1.0;
    KmCaMK = 0.15;
    aCaMK = 0.05;
    bCaMK = 0.00068;
    CaMKo = 0.05;
    KmCaM = 0.0015;
    cmdnmax_b = 0.05;
    kmcmdn = 0.00238;
    trpnmax = 0.07;
    kmtrpn = 0.0005;
    BSRmax = 0.047;
    KmBSR = 0.00087;
    BSLmax = 1.124;
    KmBSL = 0.0087;
    csqnmax = 10;
    kmcsqn = 0.8;
    cli = 24.0;
    PKNa = 0.01833;
    gkatp = 4.3195;
    fkatp = 0.0;
    K_o_n = 5;
    A_atp = 2;
    K_atp = 0.25;
    GNa = 11.7802;
    thL = 200;
    GNaL_b = 0.0279;
    Gto_b = 0.16;
    EKshift = 0;
    Kmn = 0.002;
    k2n = 500;
    PCa_b = 8.3757e-05;
    Aff = 0.6;
    tjca = 75;
    vShift = 0;
    offset = 0;
    dielConstant = 74;
    ICaL_fractionSS = 0.8;
    GKr_b = 0.0321;
    alpha_1 = 0.154375;
    beta_1 = 0.1911;
    GKs_b = 0.0011;
    GK1_b = 0.6992;
    INaCa_fractionSS = 0.35;
    kna1 = 15;
    kna2 = 5;
    kna3 = 88.12;
    kasymm = 12.5;
    wna = 6e4;
    wca = 6e4;
    wnaca = 5e3;
    kcaon = 1.5e6;
    kcaoff = 5e3;
    qna = 0.5224;
    qca = 0.167;
    KmCaAct = 150e-6;
    Gncx_b = 0.0034;
    k1p = 949.5;
    k1m = 182.4;
    k2p = 687.2;
    k2m = 39.4;
    k3p = 1899;
    k3m = 79300;
    k4p = 639;
    k4m = 40;
    Knai0 = 9.073;
    Knao0 = 27.78;
    delta = -0.155;
    Kki = 0.5;
    Kko = 0.3582;
    MgADP = 0.05;
    MgATP = 9.8;
    Kmgatp = 1.698e-7;
    H = 1e-7;
    eP = 4.2;
    Khp = 1.698e-7;
    Knap = 224;
    Kxkur = 292;
    Pnak_b = 15.4509;
    GKb_b = 0.0189;
    PNab = 1.9239e-09;
    PCab = 5.9194e-08;
    GpCa = 5e-04;
    KmCap = 0.0005;
    GClCa = 0.2843;
    GClb = 1.98e-3;
    KdClCa = 0.1;
    Fjunc = 1;
    tauNa = 2.0;
    tauK = 2.0;
    tauCa = 0.2;
    bt = 4.75;
    cajsr_half = 1.7;
    Jrel_b = 1.5378;
    Jup_b = 1.0;
    vcell = 1000.00 * 3.14000 * rad * rad * L;
    cmdnmax = (celltype == 1.00000 ? cmdnmax_b * 1.30000 : cmdnmax_b);
    ECl = ((R * T) / (zcl * F)) * log(clo / cli);
    akik = pow(ko / K_o_n, 0.240000);
    bkik = 1.00000 / (1.00000 + pow(A_atp / K_atp, 2.00000));
    thLp = 3.00000 * thL;
    GNaL = (celltype == 1.00000 ? GNaL_b * 0.600000 : GNaL_b);
    Gto = (celltype == 1.00000 ? Gto_b * 2.00000 : celltype == 2.00000 ? Gto_b * 2.00000 : Gto_b);
    Afs = 1.00000 - Aff;
    PCa = (celltype == 1.00000 ? PCa_b * 1.20000 : celltype == 2.00000 ? PCa_b * 2.00000 : PCa_b);
    Io = (0.500000 * (nao + ko + clo + 4.00000 * cao)) / 1000.00;
    GKr = (celltype == 1.00000 ? GKr_b * 1.30000 : celltype == 2.00000 ? GKr_b * 0.800000 : GKr_b);
    GKs = (celltype == 1.00000 ? GKs_b * 1.40000 : GKs_b);
    GK1 = (celltype == 1.00000 ? GK1_b * 1.20000 : celltype == 2.00000 ? GK1_b * 1.30000 : GK1_b);
    GKb = (celltype == 1.00000 ? GKb_b * 0.600000 : GKb_b);
    a_rel = (0.500000 * bt) / 1.00000;
    btp = 1.25000 * bt;
    upScale = (celltype == 1.00000 ? 1.30000 : 1.00000);
    Ageo = 2.00000 * 3.14000 * rad * rad + 2.00000 * 3.14000 * rad * L;
    PCap = 1.10000 * PCa;
    PCaNa = 0.00125000 * PCa;
    PCaK = 0.000357400 * PCa;
    constA = 1.82000e+06 * pow(dielConstant * T, -1.50000);
    a_relp = (0.500000 * btp) / 1.00000;
    Acap = 2.00000 * Ageo;
    PCaNap = 0.00125000 * PCap;
    PCaKp = 0.000357400 * PCap;
    gamma_cao = exp(-constA * 4.00000 * (pow(Io, 1.0 / 2) / (1.00000 + pow(Io, 1.0 / 2)) - 0.300000 * Io));
    gamma_nao = exp(-constA * 1.00000 * (pow(Io, 1.0 / 2) / (1.00000 + pow(Io, 1.0 / 2)) - 0.300000 * Io));
    gamma_ko = exp(-constA * 1.00000 * (pow(Io, 1.0 / 2) / (1.00000 + pow(Io, 1.0 / 2)) - 0.300000 * Io));
    vmyo = 0.680000 * vcell;
    vnsr = 0.0552000 * vcell;
    vjsr = 0.00480000 * vcell;
    vss = 0.0200000 * vcell;
    h10_i = kasymm + 1.00000 + (nao / kna1) * (1.00000 + nao / kna2);
    h11_i = (nao * nao) / (h10_i * kna1 * kna2);
    h12_i = 1.00000 / h10_i;
    k1_i = h12_i * cao * kcaon;
    k2_i = kcaoff;
    k5_i = kcaoff;
    Gncx = (celltype == 1.00000 ? Gncx_b * 1.10000 : celltype == 2.00000 ? Gncx_b * 1.40000 : Gncx_b);
    h10_ss = kasymm + 1.00000 + (nao / kna1) * (1.00000 + nao / kna2);
    h11_ss = (nao * nao) / (h10_ss * kna1 * kna2);
    h12_ss = 1.00000 / h10_ss;
    k1_ss = h12_ss * cao * kcaon;
    k2_ss = kcaoff;
    k5_ss = kcaoff;
    b1 = k1m * MgADP;
    a2 = k2p;
    a4 = ((k4p * MgATP) / Kmgatp) / (1.00000 + MgATP / Kmgatp);
    Pnak = (celltype == 1.00000 ? Pnak_b * 0.900000 : celltype == 2.00000 ? Pnak_b * 0.700000 : Pnak_b);
}
void ToRORd::set_default_initial_state(double* Y_old_)
{
    v_old_ = -88.7638;
    CaMKt_old_ = 0.0111;
    cass_old_ = 7.0305e-5;
    nai_old_ = 12.1025;
    nass_old_ = 12.1029;
    ki_old_ = 142.3002;
    kss_old_ = 142.3002;
    cansr_old_ = 1.5211;
    cajsr_old_ = 1.5214;
    cai_old_ = 8.1583e-05;
    m_old_ = 8.0572e-4;
    h_old_ = 0.8286;
    j_old_ = 0.8284;
    hp_old_ = 0.6707;
    jp_old_ = 0.8281;
    mL_old_ = 1.629e-4;
    hL_old_ = 0.5255;
    hLp_old_ = 0.2872;
    a_old_ = 9.5098e-4;
    iF_old_ = 0.9996;
    iS_old_ = 0.5936;
    ap_old_ = 4.8454e-4;
    iFp_old_ = 0.9996;
    iSp_old_ = 0.6538;
    d_old_ = 8.1084e-9;
    ff_old_ = 1.0;
    fs_old_ = 0.939;
    fcaf_old_ = 1.0;
    fcas_old_ = 0.9999;
    jca_old_ = 1.0;
    ffp_old_ = 1.0;
    fcafp_old_ = 1.0;
    nca_ss_old_ = 6.6462e-4;
    nca_i_old_ = 0.0012;
    C1_old_ = 7.0344e-4;
    C2_old_ = 8.5109e-4;
    C3_old_ = 0.9981;
    I_old_ = 1.3289e-5;
    O_old_ = 3.7585e-4;
    xs1_old_ = 0.248;
    xs2_old_ = 1.7707e-4;
    Jrel_np_old_ = 1.6129e-22;
    Jrel_p_old_ = 1.2475e-20;
}

double ToRORd::calc_stimulus(double* pars, double t)
{

}

void ToRORd::calc_algs_nl(double* algs, double* pars, double* Y_old_, double time)
{
    hLss = 1.00000 / (1.00000 + exp((v_old_ + 87.6100) / 7.48800));
    hLssp = 1.00000 / (1.00000 + exp((v_old_ + 93.8100) / 7.48800));
    jcass = 1.00000 / (1.00000 + exp((v_old_ + 18.0800) / 2.79160));
    mss = 1.00000 / pow(1.00000 + exp(-(v_old_ + 56.8600) / 9.03000), 2.00000);
    tm = 0.129200 * exp(-pow((v_old_ + 45.7900) / 15.5400, 2.00000)) + 0.0648700 * exp(-pow((v_old_ - 4.82300) / 51.1200, 2.00000));
    mLss = 1.00000 / (1.00000 + exp(-(v_old_ + 42.8500) / 5.26400));
    tmL = 0.129200 * exp(-pow((v_old_ + 45.7900) / 15.5400, 2.00000)) + 0.0648700 * exp(-pow((v_old_ - 4.82300) / 51.1200, 2.00000));
    ass = 1.00000 / (1.00000 + exp(-((v_old_ + EKshift) - 14.3400) / 14.8200));
    ta = 1.05150 / (1.00000 / (1.20890 * (1.00000 + exp(-((v_old_ + EKshift) - 18.4099) / 29.3814))) + 3.50000 / (1.00000 + exp((v_old_ + EKshift + 100.000) / 29.3814)));
    dss = (v_old_ >= 31.4978 ? 1.00000 : 1.07630 * exp(-1.00700 * exp(-0.0829000 * v_old_)));
    td = offset + 0.600000 + 1.00000 / (exp(-0.0500000 * (v_old_ + vShift + 6.00000)) + exp(0.0900000 * (v_old_ + vShift + 14.0000)));
    fss = 1.00000 / (1.00000 + exp((v_old_ + 19.5800) / 3.69600));
    tff = 7.00000 + 1.00000 / (0.00450000 * exp(-(v_old_ + 20.0000) / 10.0000) + 0.00450000 * exp((v_old_ + 20.0000) / 10.0000));
    tfs = 1000.00 + 1.00000 / (3.50000e-05 * exp(-(v_old_ + 5.00000) / 4.00000) + 3.50000e-05 * exp((v_old_ + 5.00000) / 6.00000));
    km2n = jca_old_ * 1.00000;
    anca_ss = 1.00000 / (k2n / km2n + pow(1.00000 + Kmn / cass_old_, 4.00000));
    anca_i = 1.00000 / (k2n / km2n + pow(1.00000 + Kmn / cai_old_, 4.00000));
    xs1ss = 1.00000 / (1.00000 + exp(-(v_old_ + 11.6000) / 8.93200));
    txs1 = 817.300 + 1.00000 / (0.000232600 * exp((v_old_ + 48.2800) / 17.8000) + 0.00129200 * exp(-(v_old_ + 210.000) / 230.000));
    assp = 1.00000 / (1.00000 + exp(-((v_old_ + EKshift) - 24.3400) / 14.8200));
    fcass = fss;
    tfcaf = 7.00000 + 1.00000 / (0.0400000 * exp(-(v_old_ - 4.00000) / 7.00000) + 0.0400000 * exp((v_old_ - 4.00000) / 7.00000));
    tfcas = 100.000 + 1.00000 / (0.000120000 * exp(-v_old_ / 3.00000) + 0.000120000 * exp(v_old_ / 7.00000));
    tffp = 2.50000 * tff;
    xs2ss = xs1ss;
    txs2 = 1.00000 / (0.0100000 * exp((v_old_ - 50.0000) / 20.0000) + 0.0193000 * exp(-(v_old_ + 66.5400) / 31.0000));
    CaMKb = (CaMKo * (1.00000 - CaMKt_old_)) / (1.00000 + KmCaM / cass_old_);
    hss = 1.00000 / pow(1.00000 + exp((v_old_ + 71.5500) / 7.43000), 2.00000);
    ah = (v_old_ >= -40.0000 ? 0.00000 : 0.0570000 * exp(-(v_old_ + 80.0000) / 6.80000));
    bh = (v_old_ >= -40.0000 ? 0.770000 / (0.130000 * (1.00000 + exp(-(v_old_ + 10.6600) / 11.1000))) : 2.70000 * exp(0.0790000 * v_old_) + 310000. * exp(0.348500 * v_old_));
    th = 1.00000 / (ah + bh);
    tfcafp = 2.50000 * tfcaf;
    jss = hss;
    aj = (v_old_ >= -40.0000 ? 0.00000 : ((-25428.0 * exp(0.244400 * v_old_) - 6.94800e-06 * exp(-0.0439100 * v_old_)) * (v_old_ + 37.7800)) / (1.00000 + exp(0.311000 * (v_old_ + 79.2300))));
    bj = (v_old_ >= -40.0000 ? (0.600000 * exp(0.0570000 * v_old_)) / (1.00000 + exp(-0.100000 * (v_old_ + 32.0000))) : (0.0242400 * exp(-0.0105200 * v_old_)) / (1.00000 + exp(-0.137800 * (v_old_ + 40.1400))));
    tj = 1.00000 / (aj + bj);
    hssp = 1.00000 / pow(1.00000 + exp((v_old_ + 77.5500) / 7.43000), 2.00000);
    iss = 1.00000 / (1.00000 + exp((v_old_ + EKshift + 43.9400) / 5.71100));
    delta_epi = (celltype == 1.00000 ? 1.00000 - 0.950000 / (1.00000 + exp((v_old_ + EKshift + 70.0000) / 5.00000)) : 1.00000);
    tiF_b = 4.56200 + 1.00000 / (0.393300 * exp(-(v_old_ + EKshift + 100.000) / 100.000) + 0.0800400 * exp((v_old_ + EKshift + 50.0000) / 16.5900));
    tiF = tiF_b * delta_epi;
    vfrt = (v_old_ * F) / (R * T);
    alpha = 0.116100 * exp(0.299000 * vfrt);
    beta = 0.244200 * exp(-1.60400 * vfrt);
    tjp = 1.46000 * tj;
    tiS_b = 23.6200 + 1.00000 / (0.00141600 * exp(-(v_old_ + EKshift + 96.5200) / 59.0500) + 1.78000e-08 * exp((v_old_ + EKshift + 114.100) / 8.07900));
    tiS = tiS_b * delta_epi;
    alpha_2 = 0.0578000 * exp(0.971000 * vfrt);
    beta_2 = 0.000349000 * exp(-1.06200 * vfrt);
    alpha_i = 0.253300 * exp(0.595300 * vfrt);
    beta_i = 0.0652500 * exp(-0.820900 * vfrt);
    dti_develop = 1.35400 + 0.000100000 / (exp(((v_old_ + EKshift) - 167.400) / 15.8900) + exp(-((v_old_ + EKshift) - 12.2300) / 0.215400));
    dti_recover = 1.00000 - 0.500000 / (1.00000 + exp((v_old_ + EKshift + 70.0000) / 20.0000));
    tiFp = dti_develop * dti_recover * tiF;
    tiSp = dti_develop * dti_recover * tiS;
    alpha_C2ToI = 5.20000e-05 * exp(1.52500 * vfrt);
    beta_ItoC2 = (beta_2 * beta_i * alpha_C2ToI) / (alpha_2 * alpha_i);
    f = Aff * ff_old_ + Afs * fs_old_;
    Afcaf = 0.300000 + 0.600000 / (1.00000 + exp((v_old_ - 10.0000) / 10.0000));
    Afcas = 1.00000 - Afcaf;
    fca = Afcaf * fcaf_old_ + Afcas * fcas_old_;
    fp = Aff * ffp_old_ + Afs * fs_old_;
    fcap = Afcaf * fcafp_old_ + Afcas * fcas_old_;
    vffrt = (v_old_ * F * F) / (R * T);
    Iss = (0.500000 * (nass_old_ + kss_old_ + cli + 4.00000 * cass_old_)) / 1000.00;
    gamma_cass = exp(-constA * 4.00000 * (pow(Iss, 1.0 / 2) / (1.00000 + pow(Iss, 1.0 / 2)) - 0.300000 * Iss));
    PhiCaL_ss = (4.00000 * vffrt * (gamma_cass * cass_old_ * exp(2.00000 * vfrt) - gamma_cao * cao)) / (exp(2.00000 * vfrt) - 1.00000);
    CaMKa = CaMKb + CaMKt_old_;
    fICaLp = 1.00000 / (1.00000 + KmCaMK / CaMKa);
    ICaL_ss = ICaL_fractionSS * ((1.00000 - fICaLp) * PCa * PhiCaL_ss * d_old_ * (f * (1.00000 - nca_ss_old_) + jca_old_ * fca * nca_ss_old_) + fICaLp * PCap * PhiCaL_ss * d_old_ * (fp * (1.00000 - nca_ss_old_) + jca_old_ * fcap * nca_ss_old_));
    Jrel_inf_b = ((-a_rel * ICaL_ss) / 1.00000) / (1.00000 + pow(cajsr_half / cajsr_old_, 8.00000));
    Jrel_inf = (celltype == 2.00000 ? Jrel_inf_b * 1.70000 : Jrel_inf_b);
    tau_rel_b = bt / (1.00000 + 0.0123000 / cajsr_old_);
    tau_rel = (tau_rel_b < 0.00100000 ? 0.00100000 : tau_rel_b);
    Jrel_infp_b = ((-a_relp * ICaL_ss) / 1.00000) / (1.00000 + pow(cajsr_half / cajsr_old_, 8.00000));
    Jrel_infp = (celltype == 2.00000 ? Jrel_infp_b * 1.70000 : Jrel_infp_b);
    tau_relp_b = btp / (1.00000 + 0.0123000 / cajsr_old_);
    tau_relp = (tau_relp_b < 0.00100000 ? 0.00100000 : tau_relp_b);
    EK = ((R * T) / (zk * F)) * log(ko / ki_old_);
    AiF = 1.00000 / (1.00000 + exp(((v_old_ + EKshift) - 213.600) / 151.200));
    AiS = 1.00000 - AiF;
    i = AiF * iF_old_ + AiS * iS_old_;
    ip = AiF * iFp_old_ + AiS * iSp_old_;
    fItop = 1.00000 / (1.00000 + KmCaMK / CaMKa);
    Ito = Gto * (v_old_ - EK) * ((1.00000 - fItop) * a_old_ * i + fItop * ap_old_ * ip);
    IKr = GKr * pow((ko / 5.00000), 1.0 / 2) * O_old_ * (v_old_ - EK);
    EKs = ((R * T) / (zk * F)) * log((ko + PKNa * nao) / (ki_old_ + PKNa * nai_old_));
    KsCa = 1.00000 + 0.600000 / (1.00000 + pow(3.80000e-05 / cai_old_, 1.40000));
    IKs = GKs * KsCa * xs1_old_ * xs2_old_ * (v_old_ - EKs);
    aK1 = 4.09400 / (1.00000 + exp(0.121700 * ((v_old_ - EK) - 49.9340)));
    bK1 = (15.7200 * exp(0.0674000 * ((v_old_ - EK) - 3.25700)) + exp(0.0618000 * ((v_old_ - EK) - 594.310))) / (1.00000 + exp(-0.162900 * ((v_old_ - EK) + 14.2070)));
    K1ss = aK1 / (aK1 + bK1);
    IK1 = GK1 * pow((ko / 5.00000), 1.0 / 2) * K1ss * (v_old_ - EK);
    Knao = Knao0 * exp(((1.00000 - delta) * vfrt) / 3.00000);
    a3 = (k3p * pow(ko / Kko, 2.00000)) / ((pow(1.00000 + nao / Knao, 3.00000) + pow(1.00000 + ko / Kko, 2.00000)) - 1.00000);
    P = eP / (1.00000 + H / Khp + nai_old_ / Knap + ki_old_ / Kxkur);
    b3 = (k3m * P * H) / (1.00000 + MgATP / Kmgatp);
    Knai = Knai0 * exp((delta * vfrt) / 3.00000);
    a1 = (k1p * pow(nai_old_ / Knai, 3.00000)) / ((pow(1.00000 + nai_old_ / Knai, 3.00000) + pow(1.00000 + ki_old_ / Kki, 2.00000)) - 1.00000);
    b2 = (k2m * pow(nao / Knao, 3.00000)) / ((pow(1.00000 + nao / Knao, 3.00000) + pow(1.00000 + ko / Kko, 2.00000)) - 1.00000);
    b4 = (k4m * pow(ki_old_ / Kki, 2.00000)) / ((pow(1.00000 + nai_old_ / Knai, 3.00000) + pow(1.00000 + ki_old_ / Kki, 2.00000)) - 1.00000);
    x1 = a4 * a1 * a2 + b2 * b4 * b3 + a2 * b4 * b3 + b3 * a1 * a2;
    x2 = b2 * b1 * b4 + a1 * a2 * a3 + a3 * b1 * b4 + a2 * a3 * b4;
    x3 = a2 * a3 * a4 + b3 * b2 * b1 + b2 * b1 * a4 + a3 * a4 * b1;
    x4 = b4 * b3 * b2 + a3 * a4 * a1 + b2 * a4 * a1 + b3 * b2 * a1;
    E1 = x1 / (x1 + x2 + x3 + x4);
    E2 = x2 / (x1 + x2 + x3 + x4);
    JnakNa = 3.00000 * (E1 * a3 - E2 * b3);
    E3 = x3 / (x1 + x2 + x3 + x4);
    E4 = x4 / (x1 + x2 + x3 + x4);
    JnakK = 2.00000 * (E4 * b1 - E3 * a1);
    INaK = Pnak * (zna * JnakNa + zk * JnakK);
    xkb = 1.00000 / (1.00000 + exp(-(v_old_ - 10.8968) / 23.9871));
    IKb = GKb * xkb * (v_old_ - EK);
    I_katp = fkatp * gkatp * akik * bkik * (v_old_ - EK);
    Istim = (VOI >= i_Stim_Start && VOI <= i_Stim_End && (VOI - i_Stim_Start) - floor((VOI - i_Stim_Start) / i_Stim_Period) * i_Stim_Period <= i_Stim_PulseDuration ? i_Stim_Amplitude : 0.00000);
    Ii = (0.500000 * (nai_old_ + ki_old_ + cli + 4.00000 * cai_old_)) / 1000.00;
    gamma_ki = exp(-constA * 1.00000 * (pow(Ii, 1.0 / 2) / (1.00000 + pow(Ii, 1.0 / 2)) - 0.300000 * Ii));
    PhiCaK_i = (1.00000 * vffrt * (gamma_ki * ki_old_ * exp(1.00000 * vfrt) - gamma_ko * ko)) / (exp(1.00000 * vfrt) - 1.00000);
    ICaK_i = (1.00000 - ICaL_fractionSS) * ((1.00000 - fICaLp) * PCaK * PhiCaK_i * d_old_ * (f * (1.00000 - nca_i_old_) + jca_old_ * fca * nca_i_old_) + fICaLp * PCaKp * PhiCaK_i * d_old_ * (fp * (1.00000 - nca_i_old_) + jca_old_ * fcap * nca_i_old_));
    JdiffK = (kss_old_ - ki_old_) / tauK;
    gamma_kss = exp(-constA * 1.00000 * (pow(Iss, 1.0 / 2) / (1.00000 + pow(Iss, 1.0 / 2)) - 0.300000 * Iss));
    PhiCaK_ss = (1.00000 * vffrt * (gamma_kss * kss_old_ * exp(1.00000 * vfrt) - gamma_ko * ko)) / (exp(1.00000 * vfrt) - 1.00000);
    ICaK_ss = ICaL_fractionSS * ((1.00000 - fICaLp) * PCaK * PhiCaK_ss * d_old_ * (f * (1.00000 - nca_ss_old_) + jca_old_ * fca * nca_ss_old_) + fICaLp * PCaKp * PhiCaK_ss * d_old_ * (fp * (1.00000 - nca_ss_old_) + jca_old_ * fcap * nca_ss_old_));
    ENa = ((R * T) / (zna * F)) * log(nao / nai_old_);
    fINap = 1.00000 / (1.00000 + KmCaMK / CaMKa);
    INa = GNa * (v_old_ - ENa) * pow(m_old_, 3.00000) * ((1.00000 - fINap) * h_old_ * j_old_ + fINap * hp_old_ * jp_old_);
    fINaLp = 1.00000 / (1.00000 + KmCaMK / CaMKa);
    INaL = GNaL * (v_old_ - ENa) * mL_old_ * ((1.00000 - fINaLp) * hL_old_ + fINaLp * hLp_old_);
    allo_i = 1.00000 / (1.00000 + pow(KmCaAct / cai_old_, 2.00000));
    hna = exp(qna * vfrt);
    h7_i = 1.00000 + (nao / kna3) * (1.00000 + 1.00000 / hna);
    h8_i = nao / (kna3 * hna * h7_i);
    k3pp_i = h8_i * wnaca;
    h1_i = 1.00000 + (nai_old_ / kna3) * (1.00000 + hna);
    h2_i = (nai_old_ * hna) / (kna3 * h1_i);
    k4pp_i = h2_i * wnaca;
    h4_i = 1.00000 + (nai_old_ / kna1) * (1.00000 + nai_old_ / kna2);
    h5_i = (nai_old_ * nai_old_) / (h4_i * kna1 * kna2);
    k7_i = h5_i * h2_i * wna;
    k8_i = h8_i * h11_i * wna;
    h9_i = 1.00000 / h7_i;
    k3p_i = h9_i * wca;
    k3_i = k3p_i + k3pp_i;
    hca = exp(qca * vfrt);
    h3_i = 1.00000 / h1_i;
    k4p_i = (h3_i * wca) / hca;
    k4_i = k4p_i + k4pp_i;
    h6_i = 1.00000 / h4_i;
    k6_i = h6_i * cai_old_ * kcaon;
    x1_i = k2_i * k4_i * (k7_i + k6_i) + k5_i * k7_i * (k2_i + k3_i);
    x2_i = k1_i * k7_i * (k4_i + k5_i) + k4_i * k6_i * (k1_i + k8_i);
    x3_i = k1_i * k3_i * (k7_i + k6_i) + k8_i * k6_i * (k2_i + k3_i);
    x4_i = k2_i * k8_i * (k4_i + k5_i) + k3_i * k5_i * (k1_i + k8_i);
    E1_i = x1_i / (x1_i + x2_i + x3_i + x4_i);
    E2_i = x2_i / (x1_i + x2_i + x3_i + x4_i);
    E3_i = x3_i / (x1_i + x2_i + x3_i + x4_i);
    E4_i = x4_i / (x1_i + x2_i + x3_i + x4_i);
    JncxNa_i = (3.00000 * (E4_i * k7_i - E1_i * k8_i) + E3_i * k4pp_i) - E2_i * k3pp_i;
    JncxCa_i = E2_i * k2_i - E1_i * k1_i;
    INaCa_i = (1.00000 - INaCa_fractionSS) * Gncx * allo_i * (zna * JncxNa_i + zca * JncxCa_i);
    INab = (PNab * vffrt * (nai_old_ * exp(vfrt) - nao)) / (exp(vfrt) - 1.00000);
    gamma_nai = exp(-constA * 1.00000 * (pow(Ii, 1.0 / 2) / (1.00000 + pow(Ii, 1.0 / 2)) - 0.300000 * Ii));
    PhiCaNa_i = (1.00000 * vffrt * (gamma_nai * nai_old_ * exp(1.00000 * vfrt) - gamma_nao * nao)) / (exp(1.00000 * vfrt) - 1.00000);
    ICaNa_i = (1.00000 - ICaL_fractionSS) * ((1.00000 - fICaLp) * PCaNa * PhiCaNa_i * d_old_ * (f * (1.00000 - nca_i_old_) + jca_old_ * fca * nca_i_old_) + fICaLp * PCaNap * PhiCaNa_i * d_old_ * (fp * (1.00000 - nca_i_old_) + jca_old_ * fcap * nca_i_old_));
    JdiffNa = (nass_old_ - nai_old_) / tauNa;
    allo_ss = 1.00000 / (1.00000 + pow(KmCaAct / cass_old_, 2.00000));
    h7_ss = 1.00000 + (nao / kna3) * (1.00000 + 1.00000 / hna);
    h8_ss = nao / (kna3 * hna * h7_ss);
    k3pp_ss = h8_ss * wnaca;
    h1_ss = 1.00000 + (nass_old_ / kna3) * (1.00000 + hna);
    h2_ss = (nass_old_ * hna) / (kna3 * h1_ss);
    k4pp_ss = h2_ss * wnaca;
    h4_ss = 1.00000 + (nass_old_ / kna1) * (1.00000 + nass_old_ / kna2);
    h5_ss = (nass_old_ * nass_old_) / (h4_ss * kna1 * kna2);
    k7_ss = h5_ss * h2_ss * wna;
    k8_ss = h8_ss * h11_ss * wna;
    h9_ss = 1.00000 / h7_ss;
    k3p_ss = h9_ss * wca;
    k3_ss = k3p_ss + k3pp_ss;
    h3_ss = 1.00000 / h1_ss;
    k4p_ss = (h3_ss * wca) / hca;
    k4_ss = k4p_ss + k4pp_ss;
    h6_ss = 1.00000 / h4_ss;
    k6_ss = h6_ss * cass_old_ * kcaon;
    x1_ss = k2_ss * k4_ss * (k7_ss + k6_ss) + k5_ss * k7_ss * (k2_ss + k3_ss);
    x2_ss = k1_ss * k7_ss * (k4_ss + k5_ss) + k4_ss * k6_ss * (k1_ss + k8_ss);
    x3_ss = k1_ss * k3_ss * (k7_ss + k6_ss) + k8_ss * k6_ss * (k2_ss + k3_ss);
    x4_ss = k2_ss * k8_ss * (k4_ss + k5_ss) + k3_ss * k5_ss * (k1_ss + k8_ss);
    E1_ss = x1_ss / (x1_ss + x2_ss + x3_ss + x4_ss);
    E2_ss = x2_ss / (x1_ss + x2_ss + x3_ss + x4_ss);
    E3_ss = x3_ss / (x1_ss + x2_ss + x3_ss + x4_ss);
    E4_ss = x4_ss / (x1_ss + x2_ss + x3_ss + x4_ss);
    JncxNa_ss = (3.00000 * (E4_ss * k7_ss - E1_ss * k8_ss) + E3_ss * k4pp_ss) - E2_ss * k3pp_ss;
    JncxCa_ss = E2_ss * k2_ss - E1_ss * k1_ss;
    INaCa_ss = INaCa_fractionSS * Gncx * allo_ss * (zna * JncxNa_ss + zca * JncxCa_ss);
    gamma_nass = exp(-constA * 1.00000 * (pow(Iss, 1.0 / 2) / (1.00000 + pow(Iss, 1.0 / 2)) - 0.300000 * Iss));
    PhiCaNa_ss = (1.00000 * vffrt * (gamma_nass * nass_old_ * exp(1.00000 * vfrt) - gamma_nao * nao)) / (exp(1.00000 * vfrt) - 1.00000);
    ICaNa_ss = ICaL_fractionSS * ((1.00000 - fICaLp) * PCaNa * PhiCaNa_ss * d_old_ * (f * (1.00000 - nca_ss_old_) + jca_old_ * fca * nca_ss_old_) + fICaLp * PCaNap * PhiCaNa_ss * d_old_ * (fp * (1.00000 - nca_ss_old_) + jca_old_ * fcap * nca_ss_old_));
    Jdiff = (cass_old_ - cai_old_) / tauCa;
    fJrelp = 1.00000 / (1.00000 + KmCaMK / CaMKa);
    Jrel = Jrel_b * ((1.00000 - fJrelp) * Jrel_np_old_ + fJrelp * Jrel_p_old_);
    Bcass = 1.00000 / (1.00000 + (BSRmax * KmBSR) / pow(KmBSR + cass_old_, 2.00000) + (BSLmax * KmBSL) / pow(KmBSL + cass_old_, 2.00000));
    gamma_cai = exp(-constA * 4.00000 * (pow(Ii, 1.0 / 2) / (1.00000 + pow(Ii, 1.0 / 2)) - 0.300000 * Ii));
    PhiCaL_i = (4.00000 * vffrt * (gamma_cai * cai_old_ * exp(2.00000 * vfrt) - gamma_cao * cao)) / (exp(2.00000 * vfrt) - 1.00000);
    ICaL_i = (1.00000 - ICaL_fractionSS) * ((1.00000 - fICaLp) * PCa * PhiCaL_i * d_old_ * (f * (1.00000 - nca_i_old_) + jca_old_ * fca * nca_i_old_) + fICaLp * PCap * PhiCaL_i * d_old_ * (fp * (1.00000 - nca_i_old_) + jca_old_ * fcap * nca_i_old_));
    ICaL = ICaL_ss + ICaL_i;
    ICaNa = ICaNa_ss + ICaNa_i;
    ICaK = ICaK_ss + ICaK_i;
    IpCa = (GpCa * cai_old_) / (KmCap + cai_old_);
    ICab = (PCab * 4.00000 * vffrt * (gamma_cai * cai_old_ * exp(2.00000 * vfrt) - gamma_cao * cao)) / (exp(2.00000 * vfrt) - 1.00000);
    IClCa_junc = ((Fjunc * GClCa) / (1.00000 + KdClCa / cass_old_)) * (v_old_ - ECl);
    IClCa_sl = (((1.00000 - Fjunc) * GClCa) / (1.00000 + KdClCa / cai_old_)) * (v_old_ - ECl);
    IClCa = IClCa_junc + IClCa_sl;
    IClb = GClb * (v_old_ - ECl);
    Jupnp = (upScale * 0.00542500 * cai_old_) / (cai_old_ + 0.000920000);
    Jupp = (upScale * 2.75000 * 0.00542500 * cai_old_) / ((cai_old_ + 0.000920000) - 0.000170000);
    fJupp = 1.00000 / (1.00000 + KmCaMK / CaMKa);
    Jleak = (0.00488250 * cansr_old_) / 15.0000;
    Jup = Jup_b * (((1.00000 - fJupp) * Jupnp + fJupp * Jupp) - Jleak);
    Bcai = 1.00000 / (1.00000 + (cmdnmax * kmcmdn) / pow(kmcmdn + cai_old_, 2.00000) + (trpnmax * kmtrpn) / pow(kmtrpn + cai_old_, 2.00000));
    Jtr = (cansr_old_ - cajsr_old_) / 60.0000;
    Bcajsr = 1.00000 / (1.00000 + (csqnmax * kmcsqn) / pow(kmcsqn + cajsr_old_, 2.00000));
}

void ToRORd::calc_algs_hh(double* algs, double* pars, double* Y_old_, double time)
{

}

void ToRORd::calc_rhs_nl(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
    calc_algs_nl(algs, pars, Y_old_, t);

    v_f_ = -(INa + INaL + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa_i + INaCa_ss + INaK + INab + IKb + IpCa + ICab + IClCa + IClb + I_katp + Istim);
    CaMKt_f_ = aCaMK * CaMKb * (CaMKb + CaMKt_old_) - bCaMK * CaMKt_old_;
    cass_f_ = Bcass * (((-(ICaL_ss - 2.00000 * INaCa_ss) * Acap) / (2.00000 * F * vss) + (Jrel * vjsr) / vss) - Jdiff);
    nai_f_ = (-(INa + INaL + 3.00000 * INaCa_i + ICaNa_i + 3.00000 * INaK + INab) * Acap) / (F * vmyo) + (JdiffNa * vss) / vmyo;
    nass_f_ = (-(ICaNa_ss + 3.00000 * INaCa_ss) * Acap) / (F * vss) - JdiffNa;
    ki_f_ = (-(((Ito + IKr + IKs + IK1 + IKb + I_katp + Istim) - 2.00000 * INaK) + ICaK_i) * Acap) / (F * vmyo) + (JdiffK * vss) / vmyo;
    kss_f_ = (-ICaK_ss * Acap) / (F * vss) - JdiffK;
    cansr_f_ = Jup - (Jtr * vjsr) / vnsr;
    cajsr_f_ = Bcajsr * (Jtr - Jrel);
    cai_f_ = Bcai * (((-((ICaL_i + IpCa + ICab) - 2.00000 * INaCa_i) * Acap) / (2.00000 * F * vmyo) - (Jup * vnsr) / vmyo) + (Jdiff * vss) / vmyo);
    m_f_ = (mss - m_old_) / tm;
    h_f_ = (hss - h_old_) / th;
    j_f_ = (jss - j_old_) / tj;
    hp_f_ = (hssp - hp_old_) / th;
    jp_f_ = (jss - jp_old_) / tjp;
    mL_f_ = (mLss - mL_old_) / tmL;
    hL_f_ = (hLss - hL_old_) / thL;
    hLp_f_ = (hLssp - hLp_old_) / thLp;
    a_f_ = (ass - a_old_) / ta;
    iF_f_ = (iss - iF_old_) / tiF;
    iS_f_ = (iss - iS_old_) / tiS;
    ap_f_ = (assp - ap_old_) / ta;
    iFp_f_ = (iss - iFp_old_) / tiFp;
    iSp_f_ = (iss - iSp_old_) / tiSp;
    d_f_ = (dss - d_old_) / td;
    ff_f_ = (fss - ff_old_) / tff;
    fs_f_ = (fss - fs_old_) / tfs;
    fcaf_f_ = (fcass - fcaf_old_) / tfcaf;
    fcas_f_ = (fcass - fcas_old_) / tfcas;
    jca_f_ = (jcass - jca_old_) / tjca;
    ffp_f_ = (fss - ffp_old_) / tffp;
    fcafp_f_ = (fcass - fcafp_old_) / tfcafp;
    nca_ss_f_ = anca_ss * k2n - nca_ss_old_ * km2n;
    nca_i_f_ = anca_i * k2n - nca_i_old_ * km2n;
    C1_f_ = (alpha_1 * C2_old_ + beta_2 * O_old_ + beta_ItoC2 * I_old_) - (beta_1 + alpha_2 + alpha_C2ToI) * C1_old_;
    C2_f_ = (alpha * C3_old_ + beta_1 * C1_old_) - (beta + alpha_1) * C2_old_;
    C3_f_ = beta * C2_old_ - alpha * C3_old_;
    I_f_ = (alpha_C2ToI * C1_old_ + alpha_i * O_old_) - (beta_ItoC2 + beta_i) * I_old_;
    O_f_ = (alpha_2 * C1_old_ + beta_i * I_old_) - (beta_2 + alpha_i) * O_old_;
    xs1_f_ = (xs1ss - xs1_old_) / txs1;
    xs2_f_ = (xs2ss - xs2_old_) / txs2;
    Jrel_np_f_ = (Jrel_inf - Jrel_np_old_) / tau_rel;
    Jrel_p_f_ = (Jrel_infp - Jrel_p_old_) / tau_relp;
}

void ToRORd::calc_rhs_hh(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{

}

void ToRORd::calc_rhs_mk(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{

}

void ToRORd::calc_hh_coeff(double* a, double* b, double* pars, double* algs, double* Y_old_, double t)
{

}

void ToRORd::prep_mk_transitions(double* algs, double* pars, double* Y_old_, double t)
{

}

void ToRORd::calc_mk_transitions(double** T, int mk_index, double* pars, double* algs, double* Y_old_, double t)
{

}

