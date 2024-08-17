#ifndef CONST_H
#define CONST_H

#include <TString.h>


#include "ErrorLogs.h"
#include "MainMenu.h"
#include "../../Include/prod2root_analysis_code/fort_common.h" // Linking of FORTRAN common block interfcommon
#include "../../Include/prod2root_analysis_code/fort_func.h"   // Linking of FORTRAN klspm00_lib library


  //Constants used in the analysis
  //Basic quantities
  const double cVel = 29.9792458;      // cm/ns
  const double hBar = 6.582119569E-34; // MeV*s
  const double eleCh  = 1.602176634E-19;     // C

  //Particles' masses
  const double mPhi = 1019.461;        // MeV/c^2
  const double mK0  = 497.611;          // MeV/c^2
  const double mPi0 = 134.9768;        // MeV/c^2
  const double mPiCh = 139.57039;      // MeV/c^2
  const double mMuon = 105.6583755;      // MeV/c^2
  const double mElec = 0.510998950;     // MeV/c^2

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

  const unsigned int channNum = 6;
  
  const TString channName[channNum] = {"K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}",
                                      "Regeneration",
                                      "#omega#pi^{0}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}",
                                      "K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}3#pi^{0}",
                                      "K_{S}K_{L}#rightarrow#pi^{#pm}l^{#mp}#nu#pi^{0}#pi^{0}",
                                      "Other bcg"};
  const TString dataName = "DATA";
  const TString mcSumName = "MC sum";

  const Color_t channColor[channNum] = {kRed, kGreen, kViolet, kCyan, kBlue, kGreen-1};
  const Color_t dataColor = kBlack;
  const Color_t mcSumColor = kOrange;

  TString path_tmp = "/internal/big_one/4/users/gamrat/old_root_files";

  const int
          firstFile = 1,
          lastFile = 56;

  ErrorHandling::ErrorLogs logger;
  bool 
    firstFileRangeErr, 
    lastFileRangeErr, 
    dataTypeErr,
    menuRangeErr;


  Controls::Menu mainMenu(0);

#endif
