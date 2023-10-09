#ifndef CONST_H
#define CONST_H

#include <TString.h>
#include <TLegend.h>


  //Constants used in the analysis
  //Basic quantities
  const double c_vel = 29.9792458;      // cm/ns
  const double h_bar = 6.582119569E-34; // MeV*s
  const double e_ch  = 1.602176634E-19;     // C

  //Particles' masses
  const double m_phi = 1019.461;        // MeV/c^2
  const double m_k0  = 497.611;          // MeV/c^2
  const double m_pi0 = 134.9768;        // MeV/c^2
  const double m_pich = 139.57039;      // MeV/c^2
  const double m_mu = 105.6583755;      // MeV/c^2
  const double m_ele = 0.510998950;     // MeV/c^2

  //Branching ratios
  //Phi
  const double br_phi_kskl = 0.339;
  const double br_phi_omegapi0 = 4.7E-5;

  //K-short
  const double br_ks_pi0pi0 = 0.3069;
  const double br_ks_pippim = 0.6920;
  const double br_ks_pippimgamma = 1.79E-3;
  const double br_ks_piele = 7.04E-4;
  const double br_ks_pimu = 4.56E-4;

  //K-long
  const double br_kl_pi0pi0 = 8.64E-4;
  const double br_kl_pippim = 1.967E-3;
  const double br_kl_pippimpi0 = 0.1254;
  const double br_kl_3pi0 = 0.1952;
  const double br_kl_piele = 0.4055;
  const double br_kl_pimu = 0.2704;

  //Kaons' properties and CPV
  const double tau_S_nonCPT = 0.89564E-1; // ns
  const double tau_S_CPT = 0.8954E-1;     // ns 
  const double tau_L = 51.16;           // ns
  const double delta_mass_nonCPT = 0.5289E10; // hbar s^-1
  const double delta_mass_CPT = 0.5293E10; // hbar s^-1 
  const double mod_epsilon = 2.228E-3;
  const double Re = 1.66E-3;
  const double Im_nonCPT = -0.11;   // deg
  const double Im_CPT = -0.002;     // deg
  const double phi_pm_nonCPT = 43.4;    // deg
  const double phi_pm_CPT = 43.51;      // deg
  const double phi_00_nonCPT = 43.7;    // deg
  const double phi_00_CPT = 43.52;      // deg

  //General
  const unsigned int MaxNumtrkv = 200;
  const unsigned int MaxNumVtx = 200;
  const unsigned int MIN_CLU_ENE = 20;

  const unsigned int chann_num = 6;
  
  const TString chann_name[chann_num] = {"K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}",
                                      "Regeneration",
                                      "#omega#pi^{0}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}",
                                      "K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}3#pi^{0}",
                                      "K_{S}K_{L}#rightarrow#pi^{#pm}l^{#mp}#nu#pi^{0}#pi^{0}",
                                      "Other bcg"};
  const TString data_name = "DATA";
  const TString mcsum_name = "MC sum";

  const Color_t chann_color[chann_num] = {kRed, kGreen, kViolet, kCyan, kBlue, kGreen-1};
  const Color_t data_color = kBlack;
  const Color_t mcsum_color = kOrange;

  //Paths



#endif
