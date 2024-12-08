#ifndef CONST_H
#define CONST_H

#include <json.hpp>
#include <fstream>
#include <ctime>
#include <stdlib.h>

#include <TString.h>
#include <TStyle.h>
#include <TChain.h>
#include <TLorentzVector.h>

#include <ErrorLogs.h>
#include <MainMenu.h>
#include <RealTimeIntegration.h>

using json = nlohmann::json;

// Get the env variable for properties

const std::string
      kloedataPath = getenv("KLOE_DBV26_DK0"),
      workdirPath = getenv("WORKDIR"),
      chainFiles = kloedataPath + "/*.root",
      pdgConstFilePath = workdirPath + "/scripts/Scripts/Subanalysis/Properties/pdg_api/pdg_const.json",
      propertiesPath = getenv("PROPERTIESKLOE"),
      propName = propertiesPath + "/properties.json";

static std::ifstream propertyFile(propName.c_str());
static json properties = json::parse(propertyFile);

inline const std::string currentDateTime()
{
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  return buf;
}

struct PDGids
{
  const TString 
          Re = "/S013EPS", 
          Im = "/S013EPI", 
          K0mass = "/S011M", 
          TauS = "/S012T", 
          TauL = "/S013T",
          deltaM = "/S013D",
          modEps = "/S013EP",
          phiPM = "/S013F+-",
          phi00 = "/S013FOO";

};

// const TString constName = "/internal/big_one/4/users/gamrat/scripts/Scripts/Properties/pdg_api/pdg_const.json";
// static std::ifstream fconst(constName);
// static json constants = json::parse(fconst);

// Constants used in the analysis
// Basic quantities
const double cVel = 29.9792458;       // cm/ns
const double hBar = 6.582119569E-34;  // MeV*s
const double eleCh = 1.602176634E-19; // C

// Particles' masses
const double mPhi = 1019.461;     // MeV/c^2
const double mK0 = 1;//(Double_t)constants["values"]["/S011M"];       // MeV/c^2
const double mPi0 = 134.9768;     // MeV/c^2
const double mPiCh = 139.57039;   // MeV/c^2
const double mMuon = 105.6583755; // MeV/c^2
const double mElec = 0.510998950; // MeV/c^2
const double mOmega = 782.66;     // MeV/c^2

// Branching ratios
// Phi
const double br_phi_kskl = 0.339;
const double br_phi_omegapi0 = 4.7E-5;

// K-short
const double br_ks_pi0pi0 = 0.3069;
const double br_ks_pippim = 0.6920;
const double br_ks_pippimgamma = 1.79E-3;
const double br_ks_piele = 7.04E-4;
const double br_ks_pimu = 4.56E-4;

// K-long
const double br_kl_pi0pi0 = 8.64E-4;
const double br_kl_pippim = 1.967E-3;
const double br_kl_pippimpi0 = 0.1254;
const double br_kl_3pi0 = 0.1952;
const double br_kl_piele = 0.4055;
const double br_kl_pimu = 0.2704;

// Kaons' properties and CPV
const double tau_S_nonCPT = 1;//(Double_t)constants["values"]["/S012T"] * 1E9;     // ns
const double tau_S_CPT = 1;//0.8954E-1;         // ns
const double tau_L = 1;//(Double_t)constants["values"]["/S013T"] * 1E9;                 // ns
const double delta_mass_nonCPT = 1;//(Double_t)constants["values"]["/S013D"]; // hbar s^-1
const double delta_mass_CPT = 1;//0.5293E10;    // hbar s^-1
const double mod_epsilon = 1;//(Double_t)constants["values"]["/S013EP"];
const double Re = 1;//constants["values"]["/S013EPS"];
const double Im_nonCPT = 1;//(Double_t)constants["values"]["/S013EPI"];    // deg
const double Im_CPT = 1;//-0.002;      // deg
const double phi_pm_nonCPT = 1;//(Double_t)constants["values"]["/S013F+-"]; // deg
const double phi_pm_CPT = 1;//43.51;   // deg
const double phi_00_nonCPT = 1;//(Double_t)constants["values"]["/S013FOO"]; // deg
const double phi_00_CPT = 1;//43.52;   // deg

// General
const int T0 = 2.715; // ns
const unsigned int MaxNumtrkv = 200;
const unsigned int MIN_CLU_ENE = 20;

const unsigned int channNum = 6;

const TString
    channName[channNum] = {"K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}",
                           "Regeneration",
                           "#omega#pi^{0}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}",
                           "K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}3#pi^{0}",
                           "K_{S}K_{L}#rightarrow#pi^{#pm}l^{#mp}#nu#pi^{0}#pi^{0}",
                           "Other bcg"},
    channelInt[channNum] = {"1", "3", "4", "5", "6", "7"};
const TString dataName = "DATA";
const TString mcSumName = "MC sum";

const Color_t channColor[channNum] = {kRed, kGreen, kViolet, kCyan, kBlue, kGreen - 1};
const Color_t dataColor = kBlack;
const Color_t mcSumColor = kOrange;

const TString
    base_path = "/internal/big_one/4/users/gamrat/scripts/Scripts/",
    path_tmp = "/internal/big_one/4/users/gamrat/old_root_files",
    prod2root_path_v26 = "/data/k2/DBV-26/DK0",
    ext_root = ".root",
    ext_img = ".svg",
    ext_csv = ".csv";

const TString
    gen_vars_filename = "gen_vars_",
    mctruth_filename = "mctruth_",
    neu_tri_filename = "neuvtx_tri_rec_",
    neu_triangle_filename = "neuvtx_triangle_rec_",
    neu_trilateration_kin_fit_filename = "neuvtx_tri_kin_fit_",
    omega_rec_filename = "omega_rec_",
    cut_vars_filename = "cut_vars_";

const TString
    gen_vars_dir = base_path + "GeneratedVars/",
    neutrec_dir = base_path + "Neutrec/",
    cpfit_dir = base_path + "CPFit/",
    omegarec_dir = base_path + "OmegaRec/",
    efficiency_dir = base_path + "EfficiencyAnalysis/",
    root_files_dir = "root_files/",
    input_dir = "input/",
    logs_dir = "log/",
    result_dir = "results/",
    img_dir = "img/";

const TString
    gen_vars_tree = "h_gen_vars",
    neutrec_triangle_tree = "h_triangle",
    neutrec_tri_tree = "h_tri",
    neutrec_kin_fit_tree = "h_tri_kin_fit",
    omegarec_tree = "h_omega_rec",
    cut_vars_tree = "h_cut_vars";

const int
    firstFileMax = 1,
    lastFileMax = 56;

struct NeuPart
{
  TLorentzVector
      vtxLAB,
      momentumLAB,
      vtxPhiCM,
      momentumPhiCM,
      vtxMotherCM,
      momentumMotherCM;
  Double_t
      InvMass,
      TotMomentum;

  Int_t
      CluNum[2];
};

struct ChPart
{
  TLorentzVector
      vtxLAB,
      momentumLAB,
      vtxPhiCM,
      momentumPhiCM,
      vtxMotherCM,
      momentumMotherCM;
  Double_t
      InvMass,
      TotMomentum;

  Int_t
      TrkNum,
      VtxNum;
};

struct Phi
{
  TLorentzVector
      vtxLAB,
      momentumLAB,
      vtxPhiCM,
      momentumPhiCM;

  Double_t
      InvMass,
      TotMomentum;

  Int_t
      TrkNum,
      VtxNum;
};

struct BaseKinematics
{
  Float_t
      Kchboost[9],
      Knereclor[9],
      Knerec[9],
      Kchrec[9],
      Kchmc[9],
      Knemc[9],
      ip[3],
      ipmc[3],
      phi_mom[4],
      Dtmc,
      Dtrec,
      Dtboostlor,
      Tcl[50],
      cluster[5][200],
      bhabha_vtx[3],
      T0step1,
      Chi2,
      minv4gam,
      Qmiss,
      trk[2][4],
      Pgamrec[4][4],
      omega[9],
      pi0[2][6];

  UChar_t
      mctruth,
      mcisr,
      g4taken[4],
      mcflag,
      pidmc[MaxNumtrkv],
      vtxmc[MaxNumtrkv],
      mother[MaxNumtrkv],
      Vtx[MaxNumtrkv],
      ncll[200],
      vtaken[3],
      VtxCh[3][MaxNumtrkv];

  Int_t
      ntmc,
      nvtxmc,
      nclu,
      nv,
      ntv,
      mctruth_int,
      errFlag;
};

struct NeutRec4
{
  Float_t
      Knerec[10],
      Photons[4][9],
      chi2min;
  Int_t
      gtaken[4],
      done;
};

inline void setGlobalStyle()
{
  // Global Style of histograms, pads, etc.

  gStyle->SetOptStat("iouMn");

  gStyle->SetFitFormat("6.2g");
  gStyle->SetStatFormat("6.2g");

  gStyle->SetCanvasDefH(750);
  gStyle->SetCanvasDefW(750);

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetOptLogz(1);
  gStyle->SetPalette(1);

  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(2);

  gStyle->SetLabelSize(0.04, "X");
  gStyle->SetLabelSize(0.04, "Y");
  gStyle->SetLabelSize(0.04, "Z");
  gStyle->SetLabelOffset(0.02, "X");
  gStyle->SetLabelOffset(0.02, "Y");
  gStyle->SetLabelOffset(0.02, "Z");
  gStyle->SetLabelFont(62, "X");
  gStyle->SetLabelFont(62, "Y");
  gStyle->SetLabelFont(62, "Z");

  gStyle->SetTitleSize(0.05, "X");
  gStyle->SetTitleSize(0.05, "Y");
  gStyle->SetTitleSize(0.05, "Z");
  gStyle->SetTitleOffset(1.2, "X");
  gStyle->SetTitleOffset(1.5, "Y");
  gStyle->SetTitleOffset(1.2, "Z");
  gStyle->SetTitleFont(62, "X");
  gStyle->SetTitleFont(62, "Y");
  gStyle->SetTitleFont(62, "Z");

  gStyle->SetTitle("");

  gStyle->cd();
}

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
};

inline void chain_init(TChain *chain_init, UInt_t first, UInt_t last)
{
  // char *env = "KLOEFILES";

  TString path(path_tmp);

  TString fullname = "",
          dirnamemc = "MONTE_CARLO", dirnamedata = "DATA",
          filenamemc = "mc_stream62_mccard2_", filenamedata = "data_stream42_",
          extension = ".root";

  for (Int_t i = first; i <= last; i++)
  {
    fullname = path + "/" + dirnamedata + "/" + filenamedata + i + extension;
    chain_init->Add(fullname);
  }
  for (Int_t i = first; i <= last; i++)
  {

    fullname = path + "/" + dirnamemc + "/" + filenamemc + i + extension;
    chain_init->Add(fullname);
  }
};

#endif
