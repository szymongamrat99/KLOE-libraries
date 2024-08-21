#ifndef CONST_H
#define CONST_H

#include <TString.h>

#include "Codes/ErrorLogs.h"
#include "Codes/MainMenu.h"

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
  const int T0 = 2.715; // ns
  const unsigned int MaxNumtrkv = 200;
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

  const TString
              base_path = "/internal/big_one/4/users/gamrat/scripts/Scripts/", 
              path_tmp = "/internal/big_one/4/users/gamrat/old_root_files",
              ext_root = ".root",
              ext_img = ".svg";

  const TString
              gen_vars_filename = "gen_vars_",
              neu_tri_filename = "neuvtx_tri_rec_",
              neu_triangle_filename = "neuvtx_triangle_rec_",
              neu_trilateration_kin_fit_filename = "neuvtx_tri_kin_fit_";

  const TString 
              gen_vars_dir = base_path + "GeneratedVars/",
              neutrec_dir = base_path + "Neutrec/",
              root_files_dir = "root_files/",
              logs_dir = "log/", 
              img_dir = "img/";

  const TString 
              gen_vars_tree = "h_gen_vars",
              neutrec_triangle_tree = "h_triangle",
              neutrec_tri_tree = "h_tri",
              neutrec_kin_fit_tree = "h_tri_kin_fit";

  const int
          firstFileMax = 1,
          lastFileMax = 56;

// Data flags inline function

inline void dataFlagSetter(Controls::DataType dataType, bool &dataFlag, int mcflag, int mctruth)
{
  switch (dataType)
  {
    case Controls::DataType::SIGNAL_TOT:
      dataFlag = (mcflag == 1 && (mctruth == 1 || mctruth == 2));
    case Controls::DataType::SIG_BCG:
      dataFlag = (mcflag == 1 && mctruth != 0);
    case Controls::DataType::MC_DATA:
      dataFlag = (mcflag == 0 || (mcflag = 1 && mctruth != 0));
  }
}

inline TString elapsedTimeHMS(double totalSeconds)
{
  int elapsedMinutes, elapsedHours;
  double elapsedSeconds;
  TString elapsedHMS;

  elapsedHours = int(totalSeconds / 3600.);
  elapsedMinutes = int(((totalSeconds / 3600.) - elapsedHours) * 60.);
  elapsedSeconds = ((((totalSeconds / 3600.) - elapsedHours) * 60.) - elapsedMinutes) * 60.;

  elapsedHMS = std::to_string(elapsedHours) + "h " + std::to_string(elapsedMinutes) + "min " + std::to_string(elapsedSeconds) + "s";

  return elapsedHMS;
}

#endif
