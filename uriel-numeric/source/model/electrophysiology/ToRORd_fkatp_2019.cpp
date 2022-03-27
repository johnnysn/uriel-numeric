#include "pch.h"
#include "model/electrophysiology/ToRORd_fkatp_2019.h"
#include "model/electrophysiology/ToRORd_fkatp_2019_defs.h"

ToRORd_fkatp_2019::ToRORd_fkatp_2019() : CellModel(17, 26, 0, 0, 278, 111) {

}

void ToRORd_fkatp_2019::calc_algs_nl(double* algs, double* pars, double* Y_old_, double time)
{
    int celltype = ENDO;

    Istim = calc_stimulus(pars, time);
    vcell =  1000.00*3.14000*rad*rad*L;
    cmdnmax = (celltype==EPI ?  cmdnmax_b*1.30000 : cmdnmax_b);
    ECl =  (( R*T)/( zcl*F))*log(clo/cli);
    akik = pow(ko/K_o_n, 0.240000);
    bkik = 1.00000/(1.00000+pow(A_atp/K_atp, 2.00000));
    GNaL = (celltype==EPI ?  GNaL_b*0.600000 : GNaL_b);
    Gto = (celltype==EPI ?  Gto_b*2.00000 : celltype==MID ?  Gto_b*2.00000 : Gto_b);    
    GKr = (celltype==EPI ?  GKr_b*1.30000 : celltype==MID ?  GKr_b*0.800000 : GKr_b);
    GKs = (celltype==EPI ?  GKs_b*1.40000 : GKs_b);
    GK1 = (celltype==EPI ?  GK1_b*1.20000 : celltype==MID ?  GK1_b*1.30000 : GK1_b);
    GKb = (celltype==EPI ?  GKb_b*0.600000 : GKb_b);
    upScale = (celltype==EPI ? 1.30000 : 1.00000);
    Ageo =  2.00000*3.14000*rad*rad+ 2.00000*3.14000*rad*L;
    Acap =  2.00000*Ageo;
    constA =  1.82000e+06*pow( dielConstant*T, - 1.50000);
    Io = ( 0.500000*(nao+ko+clo+ 4.00000*cao))/1000.00;
    gamma_nao = exp( - constA*1.00000*( pow(Io, 1.0 / 2)/(1.00000+ pow(Io, 1.0 / 2)) -  0.300000*Io));
    gamma_ko = exp( - constA*1.00000*( pow(Io, 1.0 / 2)/(1.00000+ pow(Io, 1.0 / 2)) -  0.300000*Io));
    vmyo =  0.680000*vcell;
    vnsr =  0.0552000*vcell;
    vjsr =  0.00480000*vcell;
    vss =  0.0200000*vcell;
    h10_i = kasymm+1.00000+ (nao/kna1)*(1.00000+nao/kna2);
    h11_i = ( nao*nao)/( h10_i*kna1*kna2);
    h12_i = 1.00000/h10_i;
    k1_i =  h12_i*cao*kcaon;
    k2_i = kcaoff;
    k5_i = kcaoff;
    Gncx = (celltype==EPI ?  Gncx_b*1.10000 : celltype==MID ?  Gncx_b*1.40000 : Gncx_b);
    h10_ss = kasymm+1.00000+ (nao/kna1)*(1.00000+nao/kna2);
    h11_ss = ( nao*nao)/( h10_ss*kna1*kna2);
    h12_ss = 1.00000/h10_ss;
    k1_ss =  h12_ss*cao*kcaon;
    k2_ss = kcaoff;
    k5_ss = kcaoff;
    b1 =  k1m*MgADP;
    a2 = k2p;
    a4 = (( k4p*MgATP)/Kmgatp)/(1.00000+MgATP/Kmgatp);
    Pnak = (celltype==EPI ?  Pnak_b*0.900000 : celltype==MID ?  Pnak_b*0.700000 : Pnak_b);
    km2n =  jca*1.00000;
    anca = 1.00000/(k2n/km2n+pow(1.00000+Kmn/cass, 4.00000));
    anca_i = 1.00000/(k2n/km2n+pow(1.00000+Kmn/cai, 4.00000));
    EK =  (( R*T)/( zk*F))*log(ko/ki);
    AiF = 1.00000/(1.00000+exp(((v+EKshift) - 213.600)/151.200));
    AiS = 1.00000 - AiF;
    i =  AiF*iF+ AiS*iS;
    ip =  AiF*iFp+ AiS*iSp;
    CaMKb = ( CaMKo*(1.00000 - CaMKt))/(1.00000+KmCaM/cass);
    CaMKa = CaMKb+CaMKt;
    fICaLp = 1.00000/(1.00000+KmCaMK/CaMKa);
    fItop = 1.00000/(1.00000+KmCaMK/CaMKa);
    Ito =  Gto*(v - EK)*( (1.00000 - fItop)*a_old*i+ fItop*ap*ip);
    IKr =  GKr* pow((ko/5.00000), 1.0 / 2)*O*(v - EK);
    EKs =  (( R*T)/( zk*F))*log((ko+ PKNa*nao)/(ki+ PKNa*nai));
    KsCa = 1.00000+0.600000/(1.00000+pow(3.80000e-05/cai, 1.40000));
    IKs =  GKs*KsCa*xs1*xs2*(v - EKs);
    aK1 = 4.09400/(1.00000+exp( 0.121700*((v - EK) - 49.9340)));
    bK1 = ( 15.7200*exp( 0.0674000*((v - EK) - 3.25700))+exp( 0.0618000*((v - EK) - 594.310)))/(1.00000+exp( - 0.162900*((v - EK)+14.2070)));
    K1ss = aK1/(aK1+bK1);
    IK1 =  GK1* pow((ko/5.00000), 1.0 / 2)*K1ss*(v - EK);
    vfrt = ( v*F)/( R*T);
    Knao =  Knao0*exp(( (1.00000 - delta)*vfrt)/3.00000);
    a3 = ( k3p*pow(ko/Kko, 2.00000))/((pow(1.00000+nao/Knao, 3.00000)+pow(1.00000+ko/Kko, 2.00000)) - 1.00000);
    P = eP/(1.00000+H/Khp+nai/Knap+ki/Kxkur);
    b3 = ( k3m*P*H)/(1.00000+MgATP/Kmgatp);
    Knai =  Knai0*exp(( delta*vfrt)/3.00000);
    a1 = ( k1p*pow(nai/Knai, 3.00000))/((pow(1.00000+nai/Knai, 3.00000)+pow(1.00000+ki/Kki, 2.00000)) - 1.00000);
    b2 = ( k2m*pow(nao/Knao, 3.00000))/((pow(1.00000+nao/Knao, 3.00000)+pow(1.00000+ko/Kko, 2.00000)) - 1.00000);
    b4 = ( k4m*pow(ki/Kki, 2.00000))/((pow(1.00000+nai/Knai, 3.00000)+pow(1.00000+ki/Kki, 2.00000)) - 1.00000);
    x1 =  a4*a1*a2+ b2*b4*b3+ a2*b4*b3+ b3*a1*a2;
    x2 =  b2*b1*b4+ a1*a2*a3+ a3*b1*b4+ a2*a3*b4;
    x3 =  a2*a3*a4+ b3*b2*b1+ b2*b1*a4+ a3*a4*b1;
    x4 =  b4*b3*b2+ a3*a4*a1+ b2*a4*a1+ b3*b2*a1;
    E1 = x1/(x1+x2+x3+x4);
    E2 = x2/(x1+x2+x3+x4);
    JnakNa =  3.00000*( E1*a3 -  E2*b3);
    E3 = x3/(x1+x2+x3+x4);
    E4 = x4/(x1+x2+x3+x4);
    JnakK =  2.00000*( E4*b1 -  E3*a1);
    INaK =  Pnak*( zna*JnakNa+ zk*JnakK);
    xkb = 1.00000/(1.00000+exp(- (v - 10.8968)/23.9871));
    IKb =  GKb*xkb*(v - EK);
    I_katp =  fkatp*gkatp*akik*bkik*(v - EK);
    Ii = ( 0.500000*(nai+ki+cli+ 4.00000*cai))/1000.00;
    gamma_ki = exp( - constA*1.00000*( pow(Ii, 1.0 / 2)/(1.00000+ pow(Ii, 1.0 / 2)) -  0.300000*Ii));
    vffrt = ( v*F*F)/( R*T);
    PhiCaK_i = ( 1.00000*vffrt*( gamma_ki*ki*exp( 1.00000*vfrt) -  gamma_ko*ko))/(exp( 1.00000*vfrt) - 1.00000);
    Afs = 1.00000 - Aff;
    f =  Aff*ff+ Afs*fs;
    PCa = (celltype==EPI ?  PCa_b*1.20000 : celltype==MID ?  PCa_b*2.00000 : PCa_b);
    PCap =  1.10000*PCa;
    PCaNa =  0.00125000*PCa;
    PCaK =  0.000357400*PCa;
    PCaNap =  0.00125000*PCap;
    PCaKp =  0.000357400*PCap;
    Afcaf = 0.300000+0.600000/(1.00000+exp((v - 10.0000)/10.0000));
    Afcas = 1.00000 - Afcaf;
    fca =  Afcaf*fcaf+ Afcas*fcas;
    fp =  Aff*ffp+ Afs*fs;
    fcap =  Afcaf*fcafp+ Afcas*fcas;
    ICaK_i =  (1.00000 - ICaL_fractionSS)*( (1.00000 - fICaLp)*PCaK*PhiCaK_i*d*( f*(1.00000 - nca_i)+ jca*fca*nca_i)+ fICaLp*PCaKp*PhiCaK_i*d*( fp*(1.00000 - nca_i)+ jca*fcap*nca_i));
    JdiffK = (kss - ki)/tauK;
    Iss = ( 0.500000*(nass+kss+cli+ 4.00000*cass))/1000.00;
    gamma_kss = exp( - constA*1.00000*( pow(Iss, 1.0 / 2)/(1.00000+ pow(Iss, 1.0 / 2)) -  0.300000*Iss));
    PhiCaK_ss = ( 1.00000*vffrt*( gamma_kss*kss*exp( 1.00000*vfrt) -  gamma_ko*ko))/(exp( 1.00000*vfrt) - 1.00000);
    ICaK_ss =  ICaL_fractionSS*( (1.00000 - fICaLp)*PCaK*PhiCaK_ss*d*( f*(1.00000 - nca_ss)+ jca*fca*nca_ss)+ fICaLp*PCaKp*PhiCaK_ss*d*( fp*(1.00000 - nca_ss)+ jca*fcap*nca_ss));
    ENa =  (( R*T)/( zna*F))*log(nao/nai);
    fINap = 1.00000/(1.00000+KmCaMK/CaMKa);
    INa =  GNa*(v - ENa)*pow(m, 3.00000)*( (1.00000 - fINap)*h*j+ fINap*hp*jp);
    fINaLp = 1.00000/(1.00000+KmCaMK/CaMKa);
    INaL =  GNaL*(v - ENa)*mL*( (1.00000 - fINaLp)*hL+ fINaLp*hLp);
    allo_i = 1.00000/(1.00000+pow(KmCaAct/cai, 2.00000));
    hna = exp( qna*vfrt);
    h7_i = 1.00000+ (nao/kna3)*(1.00000+1.00000/hna);
    h8_i = nao/( kna3*hna*h7_i);
    k3pp_i =  h8_i*wnaca;
    h1_i = 1.00000+ (nai/kna3)*(1.00000+hna);
    h2_i = ( nai*hna)/( kna3*h1_i);
    k4pp_i =  h2_i*wnaca;
    h4_i = 1.00000+ (nai/kna1)*(1.00000+nai/kna2);
    h5_i = ( nai*nai)/( h4_i*kna1*kna2);
    k7_i =  h5_i*h2_i*wna;
    k8_i =  h8_i*h11_i*wna;
    h9_i = 1.00000/h7_i;
    k3p_i =  h9_i*wca;
    k3_i = k3p_i+k3pp_i;
    hca = exp( qca*vfrt);
    h3_i = 1.00000/h1_i;
    k4p_i = ( h3_i*wca)/hca;
    k4_i = k4p_i+k4pp_i;
    h6_i = 1.00000/h4_i;
    k6_i =  h6_i*cai*kcaon;
    x1_i =  k2_i*k4_i*(k7_i+k6_i)+ k5_i*k7_i*(k2_i+k3_i);
    x2_i =  k1_i*k7_i*(k4_i+k5_i)+ k4_i*k6_i*(k1_i+k8_i);
    x3_i =  k1_i*k3_i*(k7_i+k6_i)+ k8_i*k6_i*(k2_i+k3_i);
    x4_i =  k2_i*k8_i*(k4_i+k5_i)+ k3_i*k5_i*(k1_i+k8_i);
    E1_i = x1_i/(x1_i+x2_i+x3_i+x4_i);
    E2_i = x2_i/(x1_i+x2_i+x3_i+x4_i);
    E3_i = x3_i/(x1_i+x2_i+x3_i+x4_i);
    E4_i = x4_i/(x1_i+x2_i+x3_i+x4_i);
    JncxNa_i = ( 3.00000*( E4_i*k7_i -  E1_i*k8_i)+ E3_i*k4pp_i) -  E2_i*k3pp_i;
    JncxCa_i =  E2_i*k2_i -  E1_i*k1_i;
    INaCa_i =  (1.00000 - INaCa_fractionSS)*Gncx*allo_i*( zna*JncxNa_i+ zca*JncxCa_i);
    INab = ( PNab*vffrt*( nai*exp(vfrt) - nao))/(exp(vfrt) - 1.00000);
    gamma_nai = exp( - constA*1.00000*( pow(Ii, 1.0 / 2)/(1.00000+ pow(Ii, 1.0 / 2)) -  0.300000*Ii));
    PhiCaNa_i = ( 1.00000*vffrt*( gamma_nai*nai*exp( 1.00000*vfrt) -  gamma_nao*nao))/(exp( 1.00000*vfrt) - 1.00000);
    ICaNa_i =  (1.00000 - ICaL_fractionSS)*( (1.00000 - fICaLp)*PCaNa*PhiCaNa_i*d*( f*(1.00000 - nca_i)+ jca*fca*nca_i)+ fICaLp*PCaNap*PhiCaNa_i*d*( fp*(1.00000 - nca_i)+ jca*fcap*nca_i));
    JdiffNa = (nass - nai)/tauNa;
    allo_ss = 1.00000/(1.00000+pow(KmCaAct/cass, 2.00000));
    h7_ss = 1.00000+ (nao/kna3)*(1.00000+1.00000/hna);
    h8_ss = nao/( kna3*hna*h7_ss);
    k3pp_ss =  h8_ss*wnaca;
    h1_ss = 1.00000+ (nass/kna3)*(1.00000+hna);
    h2_ss = ( nass*hna)/( kna3*h1_ss);
    k4pp_ss =  h2_ss*wnaca;
    h4_ss = 1.00000+ (nass/kna1)*(1.00000+nass/kna2);
    h5_ss = ( nass*nass)/( h4_ss*kna1*kna2);
    k7_ss =  h5_ss*h2_ss*wna;
    k8_ss =  h8_ss*h11_ss*wna;
    h9_ss = 1.00000/h7_ss;
    k3p_ss =  h9_ss*wca;
    k3_ss = k3p_ss+k3pp_ss;
    h3_ss = 1.00000/h1_ss;
    k4p_ss = ( h3_ss*wca)/hca;
    k4_ss = k4p_ss+k4pp_ss;
    h6_ss = 1.00000/h4_ss;
    k6_ss =  h6_ss*cass*kcaon;
    x1_ss =  k2_ss*k4_ss*(k7_ss+k6_ss)+ k5_ss*k7_ss*(k2_ss+k3_ss);
    x2_ss =  k1_ss*k7_ss*(k4_ss+k5_ss)+ k4_ss*k6_ss*(k1_ss+k8_ss);
    x3_ss =  k1_ss*k3_ss*(k7_ss+k6_ss)+ k8_ss*k6_ss*(k2_ss+k3_ss);
    x4_ss =  k2_ss*k8_ss*(k4_ss+k5_ss)+ k3_ss*k5_ss*(k1_ss+k8_ss);
    E1_ss = x1_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
    E2_ss = x2_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
    E3_ss = x3_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
    E4_ss = x4_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
    JncxNa_ss = ( 3.00000*( E4_ss*k7_ss -  E1_ss*k8_ss)+ E3_ss*k4pp_ss) -  E2_ss*k3pp_ss;
    JncxCa_ss =  E2_ss*k2_ss -  E1_ss*k1_ss;
    INaCa_ss =  INaCa_fractionSS*Gncx*allo_ss*( zna*JncxNa_ss+ zca*JncxCa_ss);
    gamma_nass = exp( - constA*1.00000*( pow(Iss, 1.0 / 2)/(1.00000+ pow(Iss, 1.0 / 2)) -  0.300000*Iss));
    PhiCaNa_ss = ( 1.00000*vffrt*( gamma_nass*nass*exp( 1.00000*vfrt) -  gamma_nao*nao))/(exp( 1.00000*vfrt) - 1.00000);
    ICaNa_ss =  ICaL_fractionSS*( (1.00000 - fICaLp)*PCaNa*PhiCaNa_ss*d*( f*(1.00000 - nca_ss)+ jca*fca*nca_ss)+ fICaLp*PCaNap*PhiCaNa_ss*d*( fp*(1.00000 - nca_ss)+ jca*fcap*nca_ss));
    Jdiff = (cass - cai)/tauCa;
    fJrelp = 1.00000/(1.00000+KmCaMK/CaMKa);
    Jrel =  Jrel_b*( (1.00000 - fJrelp)*Jrel_np+ fJrelp*Jrel_p);
    Bcass = 1.00000/(1.00000+( BSRmax*KmBSR)/pow(KmBSR+cass, 2.00000)+( BSLmax*KmBSL)/pow(KmBSL+cass, 2.00000));
    gamma_cai = exp( - constA*4.00000*( pow(Ii, 1.0 / 2)/(1.00000+ pow(Ii, 1.0 / 2)) -  0.300000*Ii));
    Io = ( 0.500000*(nao+ko+clo+ 4.00000*cao))/1000.00;
    gamma_cao = exp( - constA*4.00000*( pow(Io, 1.0 / 2)/(1.00000+ pow(Io, 1.0 / 2)) -  0.300000*Io));
    PhiCaL_i = ( 4.00000*vffrt*( gamma_cai*cai*exp( 2.00000*vfrt) -  gamma_cao*cao))/(exp( 2.00000*vfrt) - 1.00000);
    ICaL_i =  (1.00000 - ICaL_fractionSS)*( (1.00000 - fICaLp)*PCa*PhiCaL_i*d*( f*(1.00000 - nca_i)+ jca*fca*nca_i)+ fICaLp*PCap*PhiCaL_i*d*( fp*(1.00000 - nca_i)+ jca*fcap*nca_i));
    gamma_cass = exp( - constA*4.00000*( pow(Iss, 1.0 / 2)/(1.00000+ pow(Iss, 1.0 / 2)) -  0.300000*Iss));
    PhiCaL_ss = ( 4.00000*vffrt*( gamma_cass*cass*exp( 2.00000*vfrt) -  gamma_cao*cao))/(exp( 2.00000*vfrt) - 1.00000);
    ICaL_ss =  ICaL_fractionSS*( (1.00000 - fICaLp)*PCa*PhiCaL_ss*d*( f*(1.00000 - nca_ss)+ jca*fca*nca_ss)+ fICaLp*PCap*PhiCaL_ss*d*( fp*(1.00000 - nca_ss)+ jca*fcap*nca_ss));
    ICaL = ICaL_ss+ICaL_i;
    ICaNa = ICaNa_ss+ICaNa_i;
    ICaK = ICaK_ss+ICaK_i;
    IpCa = ( GpCa*cai)/(KmCap+cai);
    ICab = ( PCab*4.00000*vffrt*( gamma_cai*cai*exp( 2.00000*vfrt) -  gamma_cao*cao))/(exp( 2.00000*vfrt) - 1.00000);
    IClCa_junc =  (( Fjunc*GClCa)/(1.00000+KdClCa/cass))*(v - ECl);
    IClCa_sl =  (( (1.00000 - Fjunc)*GClCa)/(1.00000+KdClCa/cai))*(v - ECl);
    IClCa = IClCa_junc+IClCa_sl;
    IClb =  GClb*(v - ECl);
    Jupnp = ( upScale*0.00542500*cai)/(cai+0.000920000);
    Jupp = ( upScale*2.75000*0.00542500*cai)/((cai+0.000920000) - 0.000170000);
    fJupp = 1.00000/(1.00000+KmCaMK/CaMKa);
    Jleak = ( 0.00488250*cansr)/15.0000;
    Jup =  Jup_b*(( (1.00000 - fJupp)*Jupnp+ fJupp*Jupp) - Jleak);
    Bcai = 1.00000/(1.00000+( cmdnmax*kmcmdn)/pow(kmcmdn+cai, 2.00000)+( trpnmax*kmtrpn)/pow(kmtrpn+cai, 2.00000));
    Jtr = (cansr - cajsr)/60.0000;
    Bcajsr = 1.00000/(1.00000+( csqnmax*kmcsqn)/pow(kmcsqn+cajsr, 2.00000));
    alpha_2 =  0.0578000*exp( 0.971000*vfrt);
    beta_2 =  0.000349000*exp( - 1.06200*vfrt);
    alpha_i =  0.253300*exp( 0.595300*vfrt);
    beta_i =  0.0652500*exp( - 0.820900*vfrt);
    alpha_C2ToI =  5.20000e-05*exp( 1.52500*vfrt);
    beta_ItoC2 = ( beta_2*beta_i*alpha_C2ToI)/( alpha_2*alpha_i);
    alpha =  0.116100*exp( 0.299000*vfrt);
    beta =  0.244200*exp( - 1.60400*vfrt);
}

void ToRORd_fkatp_2019::calc_algs_hh(double* algs, double* pars, double* Y_old_, double time)
{
    // celltype:> ENDO = 0, MID = 1, EPI = 2
    int celltype = ENDO;

    mss = 1.00000/pow(1.00000+exp(- (v+56.8600)/9.03000), 2.00000);
    tm =  0.129200*exp(- pow((v+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((v - 4.82300)/51.1200, 2.00000));

    hss = 1.00000/pow(1.00000+exp((v+71.5500)/7.43000), 2.00000);
    hssp = 1.00000/pow(1.00000+exp((v+77.5500)/7.43000), 2.00000);
    hLss = 1.00000/(1.00000+exp((v+87.6100)/7.48800));
    hLssp = 1.00000/(1.00000+exp((v+93.8100)/7.48800));
    ah = (v>=- 40.0000 ? 0.00000 :  0.0570000*exp(- (v+80.0000)/6.80000));
    bh = (v>=- 40.0000 ? 0.770000/( 0.130000*(1.00000+exp(- (v+10.6600)/11.1000))) :  2.70000*exp( 0.0790000*v)+ 310000.*exp( 0.348500*v));
    th = 1.00000/(ah+bh);
    thLp =  3.00000*thL;
    
    jss = hss;
    aj = (v>=- 40.0000 ? 0.00000 : ( ( - 25428.0*exp( 0.244400*v) -  6.94800e-06*exp( - 0.0439100*v))*(v+37.7800))/(1.00000+exp( 0.311000*(v+79.2300))));
    bj = (v>=- 40.0000 ? ( 0.600000*exp( 0.0570000*v))/(1.00000+exp( - 0.100000*(v+32.0000))) : ( 0.0242400*exp( - 0.0105200*v))/(1.00000+exp( - 0.137800*(v+40.1400))));
    tj = 1.00000/(aj+bj);
    tjp =  1.46000*tj;
    
    mLss = 1.00000/(1.00000+exp(- (v+42.8500)/5.26400));
    tmL =  0.129200*exp(- pow((v+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((v - 4.82300)/51.1200, 2.00000));

    ass = 1.00000/(1.00000+exp(- ((v+EKshift) - 14.3400)/14.8200));
    assp = 1.00000/(1.00000+exp(- ((v+EKshift) - 24.3400)/14.8200));
    ta = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- ((v+EKshift) - 18.4099)/29.3814)))+3.50000/(1.00000+exp((v+EKshift+100.000)/29.3814)));
    
    iss = 1.00000/(1.00000+exp((v+EKshift+43.9400)/5.71100));
    delta_epi = (celltype==EPI ? 1.00000 - 0.950000/(1.00000+exp((v+EKshift+70.0000)/5.00000)) : 1.00000);
    tiF_b = 4.56200+1.00000/( 0.393300*exp(- (v+EKshift+100.000)/100.000)+ 0.0800400*exp((v+EKshift+50.0000)/16.5900));
    tiF =  tiF_b*delta_epi;
    tiS_b = 23.6200+1.00000/( 0.00141600*exp(- (v+EKshift+96.5200)/59.0500)+ 1.78000e-08*exp((v+EKshift+114.100)/8.07900));
    tiS =  tiS_b*delta_epi;

    dti_develop = 1.35400+0.000100000/(exp(((v+EKshift) - 167.400)/15.8900)+exp(- ((v+EKshift) - 12.2300)/0.215400));
    dti_recover = 1.00000 - 0.500000/(1.00000+exp((v+EKshift+70.0000)/20.0000));
    tiFp =  dti_develop*dti_recover*tiF;
    tiSp =  dti_develop*dti_recover*tiS;

    dss = (v>=31.4978 ? 1.00000 :  1.07630*exp( - 1.00700*exp( - 0.0829000*v)));
    td = offset+0.600000+1.00000/(exp( - 0.0500000*(v+vShift+6.00000))+exp( 0.0900000*(v+vShift+14.0000)));

    fss = 1.00000/(1.00000+exp((v+19.5800)/3.69600));
    tff = 7.00000+1.00000/( 0.00450000*exp(- (v+20.0000)/10.0000)+ 0.00450000*exp((v+20.0000)/10.0000));
    tfs = 1000.00+1.00000/( 3.50000e-05*exp(- (v+5.00000)/4.00000)+ 3.50000e-05*exp((v+5.00000)/6.00000));
    
    fcass = fss;
    tfcaf = 7.00000+1.00000/( 0.0400000*exp(- (v - 4.00000)/7.00000)+ 0.0400000*exp((v - 4.00000)/7.00000));
    tfcas = 100.000+1.00000/( 0.000120000*exp(- v/3.00000)+ 0.000120000*exp(v/7.00000));
    tffp =  2.50000*tff;
    tfcafp =  2.50000*tfcaf;
    
    jcass = 1.00000/(1.00000+exp((v+18.0800)/2.79160));

    xs1ss = 1.00000/(1.00000+exp(- (v+11.6000)/8.93200));
    txs1 = 817.300+1.00000/( 0.000232600*exp((v+48.2800)/17.8000)+ 0.00129200*exp(- (v+210.000)/230.000));
    
    xs2ss = xs1ss;
    txs2 = 1.00000/( 0.0100000*exp((v - 50.0000)/20.0000)+ 0.0193000*exp(- (v+66.5400)/31.0000));
    
    CaMKb = ( CaMKo*(1.00000 - CaMKt))/(1.00000+KmCaM/cass);
    CaMKa = CaMKb+CaMKt;
    fICaLp = 1.00000/(1.00000+KmCaMK/CaMKa);
    PCa = (celltype==EPI ?  PCa_b*1.20000 : celltype==MID ?  PCa_b*2.00000 : PCa_b);
    vffrt = ( v*F*F)/( R*T);
    Iss = ( 0.500000*(nass+kss+cli+ 4.00000*cass))/1000.00;
    constA =  1.82000e+06*pow( dielConstant*T, - 1.50000);
    gamma_cass = exp( - constA*4.00000*( pow(Iss, 1.0 / 2)/(1.00000+ pow(Iss, 1.0 / 2)) -  0.300000*Iss));
    vfrt = ( v*F)/( R*T);
    alpha =  0.116100*exp( 0.299000*vfrt);
    beta =  0.244200*exp( - 1.60400*vfrt);
    alpha_2 =  0.0578000*exp( 0.971000*vfrt);
    beta_2 =  0.000349000*exp( - 1.06200*vfrt);
    alpha_i =  0.253300*exp( 0.595300*vfrt);
    beta_i =  0.0652500*exp( - 0.820900*vfrt);
    alpha_C2ToI =  5.20000e-05*exp( 1.52500*vfrt);
    beta_ItoC2 = ( beta_2*beta_i*alpha_C2ToI)/( alpha_2*alpha_i);
    Io = ( 0.500000*(nao+ko+clo+ 4.00000*cao))/1000.00;
    gamma_cao = exp( - constA*4.00000*( pow(Io, 1.0 / 2)/(1.00000+ pow(Io, 1.0 / 2)) -  0.300000*Io));
    PhiCaL_ss = ( 4.00000*vffrt*( gamma_cass*cass*exp( 2.00000*vfrt) -  gamma_cao*cao))/(exp( 2.00000*vfrt) - 1.00000);
    Afs = 1.00000 - Aff;
    f =  Aff*ff+ Afs*fs;
    PCa = (celltype==EPI ?  PCa_b*1.20000 : celltype==MID ?  PCa_b*2.00000 : PCa_b);
    PCap =  1.10000*PCa;
    PCaNa =  0.00125000*PCa;
    PCaK =  0.000357400*PCa;
    PCaNap =  0.00125000*PCap;
    PCaKp =  0.000357400*PCap;
    Afcaf = 0.300000+0.600000/(1.00000+exp((v - 10.0000)/10.0000));
    Afcas = 1.00000 - Afcaf;
    fca =  Afcaf*fcaf+ Afcas*fcas;
    fp =  Aff*ffp+ Afs*fs;
    fcap =  Afcaf*fcafp+ Afcas*fcas;
    ICaL_ss =  ICaL_fractionSS*( (1.00000 - fICaLp)*PCa*PhiCaL_ss*d*( f*(1.00000 - nca_ss)+ jca*fca*nca_ss)+ fICaLp*PCap*PhiCaL_ss*d*( fp*(1.00000 - nca_ss)+ jca*fcap*nca_ss));
    a_rel = ( 0.500000*bt)/1.00000;
    Jrel_inf_b = (( - a_rel*ICaL_ss)/1.00000)/(1.00000+pow(cajsr_half/cajsr, 8.00000));
    Jrel_inf = (celltype==MID ?  Jrel_inf_b*1.70000 : Jrel_inf_b);
    
    tau_rel_b = bt/(1.00000+0.0123000/cajsr);
    tau_rel = (tau_rel_b<0.00100000 ? 0.00100000 : tau_rel_b);
    btp =  1.25000*bt;
    a_relp = ( 0.500000*btp)/1.00000;
    Jrel_infp_b = (( - a_relp*ICaL_ss)/1.00000)/(1.00000+pow(cajsr_half/cajsr, 8.00000));
    Jrel_infp = (celltype==MID ?  Jrel_infp_b*1.70000 : Jrel_infp_b);
    tau_relp_b = btp/(1.00000+0.0123000/cajsr);
    tau_relp = (tau_relp_b<0.00100000 ? 0.00100000 : tau_relp_b);
}

void ToRORd_fkatp_2019::calc_rhs_nl(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
	calc_algs_nl(algs, pars, Y_old_, t);

	V_f_ = - (INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+IClCa+IClb+I_katp+Istim);
    CaMKt_f_ =  aCaMK*CaMKb*(CaMKb+CaMKt) -  bCaMK*CaMKt;
    cass_f_ =  Bcass*((( - (ICaL_ss -  2.00000*INaCa_ss)*Acap)/( 2.00000*F*vss)+( Jrel*vjsr)/vss) - Jdiff);
    nai_f_ = ( - (INa+INaL+ 3.00000*INaCa_i+ICaNa_i+ 3.00000*INaK+INab)*Acap)/( F*vmyo)+( JdiffNa*vss)/vmyo;
    nass_f_ = ( - (ICaNa_ss+ 3.00000*INaCa_ss)*Acap)/( F*vss) - JdiffNa;
    ki_f_ = ( - (((Ito+IKr+IKs+IK1+IKb+I_katp+Istim) -  2.00000*INaK)+ICaK_i)*Acap)/( F*vmyo)+( JdiffK*vss)/vmyo;
    kss_f_ = ( - ICaK_ss*Acap)/( F*vss) - JdiffK;
    cansr_f_ = Jup - ( Jtr*vjsr)/vnsr;
    cajsr_f_ =  Bcajsr*(Jtr - Jrel);
    cai_f_ =  Bcai*((( - ((ICaL_i+IpCa+ICab) -  2.00000*INaCa_i)*Acap)/( 2.00000*F*vmyo) - ( Jup*vnsr)/vmyo)+( Jdiff*vss)/vmyo);
    nca_ss_f_ =  anca*k2n -  nca_ss*km2n;
    nca_i_f_ =  anca_i*k2n -  nca_i*km2n;
    C1_f_ = ( alpha_1*C2+ beta_2*O+ beta_ItoC2*I) -  (beta_1+alpha_2+alpha_C2ToI)*C1;
    C2_f_ = ( alpha*C3+ beta_1*C1) -  (beta+alpha_1)*C2;
    C3_f_ =  beta*C2 -  alpha*C3;
    I_f_ = ( alpha_C2ToI*C1+ alpha_i*O) -  (beta_ItoC2+beta_i)*I;
    O_f_ = ( alpha_2*C1+ beta_i*I) -  (beta_2+alpha_i)*O;
}

void ToRORd_fkatp_2019::calc_rhs_hh(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
	calc_algs_hh(algs, pars, Y_old_, t);

	m_f_ = (mss - m)/tm;
    h_f_ = (hss - h)/th;
    j_f_ = (jss - j)/tj;
    hp_f_ = (hssp - hp)/th;
    jp_f_ = (jss - jp)/tjp;
    mL_f_ = (mLss - mL)/tmL;
    hL_f_ = (hLss - hL)/thL;
    hLp_f_ = (hLssp - hLp)/thLp;
    a_f_ = (ass - a_old)/ta;
    iF_f_ = (iss - iF)/tiF;
    iS_f_ = (iss - iS)/tiS;
    ap_f_ = (assp - ap)/ta;
    iFp_f_ = (iss - iFp)/tiFp;
    iSp_f_ = (iss - iSp)/tiSp;
    d_f_ = (dss - d)/td;
    ff_f_ = (fss - ff)/tff;
    fs_f_ = (fss - fs)/tfs;
    fcaf_f_ = (fcass - fcaf)/tfcaf;
    fcas_f_ = (fcass - fcas)/tfcas;
    jca_f_ = (jcass - jca)/tjca;
    ffp_f_ = (fss - ffp)/tffp;
    fcafp_f_ = (fcass - fcafp)/tfcafp;
    xs1_f_ = (xs1ss - xs1)/txs1;
    xs2_f_ = (xs2ss - xs2)/txs2;
    Jrel_np_f_ = (Jrel_inf - Jrel_np)/tau_rel;
    Jrel_p_f_ = (Jrel_infp - Jrel_p)/tau_relp;
}

void ToRORd_fkatp_2019::calc_rhs_mk(double* rhs, double* pars, double* algs, double* Y_old_, double t)
{
}

void ToRORd_fkatp_2019::calc_hh_coeff(double* a, double* b, double* pars, double* algs, double* Y_old_, double t)
{
	calc_algs_hh(algs, pars, Y_old_, t);

    m_a_ = -1.0 / tm;
    h_a_ = -1.0 / th;
    j_a_ = -1.0 / tj;
    hp_a_ = -1.0 / th;
    jp_a_ = -1.0 / tjp;
    mL_a_ = -1.0 / tmL;
    hL_a_ = -1.0 / thL;
    hLp_a_ = -1.0 / thLp;
    a_a_ = -1.0 / ta;
    iF_a_ = -1.0 / tiF;
    iS_a_ = -1.0 / tiS;
    ap_a_ = -1.0 / ta;
    iFp_a_ = -1.0 / tiFp;
    iSp_a_ = -1.0 / tiSp;
    d_a_ = -1.0 / td;
    ff_a_ = -1.0 / tff;
    fs_a_ = -1.0 / tfs;
    fcaf_a_ = -1.0 / tfcaf;
    fcas_a_ = -1.0 / tfcas;
    jca_a_ = -1.0 / tjca;
    ffp_a_ = -1.0 / tffp;
    fcafp_a_ = -1.0 / tfcafp;
    xs1_a_ = -1.0 / txs1;
    xs2_a_ = -1.0 / txs2;
    Jrel_np_a_ = -1.0 / tau_rel;
    Jrel_p_a_ = -1.0 / tau_relp;

    m_b_ = mss / tm;
    h_b_ = hss / th;
    j_b_ = jss / tj;
    hp_b_ = hssp / th;
    jp_b_ = jss / tjp;
    mL_b_ = mLss / tmL;
    hL_b_ = hLss / thL;
    hLp_b_ = hLssp / thLp;
    a_b_ = ass / ta;
    iF_b_ = iss / tiF;
    iS_b_ = iss / tiS;
    ap_b_ = assp / ta;
    iFp_b_ = iss / tiFp;
    iSp_b_ = iss / tiSp;
    d_b_ = dss / td;
    ff_b_ = fss / tff;
    fs_b_ = fss / tfs;
    fcaf_b_ = fcass / tfcaf;
    fcas_b_ = fcass / tfcas;
    jca_b_ = jcass / tjca;
    ffp_b_ = fss / tffp;
    fcafp_b_ = fcass / tfcafp;
    xs1_b_ = xs1ss / txs1;
    xs2_b_ = xs2ss / txs2;
    Jrel_np_b_ = Jrel_inf / tau_rel;
    Jrel_p_b_ = Jrel_infp / tau_relp;
}

void ToRORd_fkatp_2019::set_default_parameters(double* pars)
{
	stim_state = 0;
	stim_amplitude = -5.30e+01;
	stim_period = 1.0e+03;
	stim_start = 0.0e+00;
	stim_duration = 1.0e+00;
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
}

void ToRORd_fkatp_2019::set_default_initial_state(double* Y_old_)
{
    // celltype = ENDO
	v = -88.7638;
    CaMKt = 0.0111;
    cass = 7.0305e-5;
    nai = 12.1025;
    nass = 12.1029;
    ki = 142.3002;
    kss = 142.3002;
    cansr = 1.5211;
    cajsr = 1.5214;
    cai = 8.1583e-05;
    m = 8.0572e-4;
    h = 0.8286;
    j = 0.8284;
    hp = 0.6707;
    jp = 0.8281;
    mL = 1.629e-4;
    hL = 0.5255;
    hLp = 0.2872;
    a_old = 9.5098e-4;
    iF = 0.9996;
    iS = 0.5936;
    ap = 4.8454e-4;
    iFp = 0.9996;
    iSp = 0.6538;
    d = 8.1084e-9;
    ff = 1.0;
    fs = 0.939;
    fcaf = 1.0;
    fcas = 0.9999;
    jca = 1.0;
    ffp = 1.0;
    fcafp = 1.0;
    nca_ss = 6.6462e-4;
    nca_i = 0.0012;
    C1 = 7.0344e-4;
    C2 = 8.5109e-4;
    C3 = 0.9981;
    I = 1.3289e-5;
    O = 3.7585e-4;
    xs1 = 0.248;
    xs2 = 1.7707e-4;
    Jrel_np = 1.6129e-22;
    Jrel_p = 1.2475e-20;
}

double ToRORd_fkatp_2019::calc_stimulus(double* pars, double t)
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

void ToRORd_fkatp_2019::prep_mk_transitions(double* algs, double* pars, double* Y_old_, double t) {}
void ToRORd_fkatp_2019::calc_mk_transitions(double** Tr, int mk_index, double* pars, double* algs, double* Y_old_, double t) {}
