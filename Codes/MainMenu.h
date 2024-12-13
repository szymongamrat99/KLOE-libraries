#ifndef MAIN_MENU_H
#define MAIN_MENU_H

#include <iostream>
#include <map>

#include <TString.h>

namespace Controls
{
  enum class MainMenu
  {
    GEN_VARS = 1,
    KCHREC_NO_BOOST = 2,
    KCHREC_BOOST = 3,
    IP_EV_BY_EV = 4,
    OMEGA_REC = 5,
    KNEREC_TRILAT = 6,
    KNEREC_TRIANGLE = 7,
    EFF_SIG_TO_BCG = 8,
    KIN_FITS = 9,
    TRANSF_TO_CM = 10,
    CPV_NORM = 11,
    EXIT = 12,

    OPT_TOT = 13
  };

  enum class NeutRecMenu
  {
    BARE_TRILATERATION = 1,
    TRILATERATION_KIN_FIT = 2,
    TRIANGLE = 3,
    COMP_OF_MET = 4,
    EXIT = 5,

    OPT_TOT = 6
  };

  enum class CPFitMenu
  {
    HALF_SIGNAL_MC = 1,
    HALF_SIG_BCG_MC = 2,
    MC_DATA = 3,
    EXIT = 4,

    OPT_TOT = 5
  };

  enum class GenVars
  {
    GEN_VARS = 1,
    SPLIT_CHANN = 2,
    EXIT = 3,

    OPT_TOT = 4
  };

  enum class OmegaRec
  {
    OMEGA_REC = 1,
    OMEGA_CUTS = 2,
    PLOTS = 3,
    EXIT = 4,

    OPT_TOT = 5
  };

  enum class Efficiency
  {
    EFF_SCAN = 1,
    EFF_DIST = 2,
    EXIT = 3,

    OPT_TOT = 4
  };

  enum class DataType
  {
    SIGNAL_TOT = 1,
    SIG_BCG = 2,
    MC_DATA = 3,

    OPT_TOT = 4
  };

  inline std::istream &operator>>(std::istream &is, MainMenu &mainMenuOpt)
  {
    int a;
    is >> a;
    mainMenuOpt = static_cast<MainMenu>(a);

    return is;
  }

  inline std::istream &operator>>(std::istream &is, NeutRecMenu &NeutRecMenuOpt)
  {
    int a;
    is >> a;
    NeutRecMenuOpt = static_cast<NeutRecMenu>(a);

    return is;
  }

  inline std::istream &operator>>(std::istream &is, CPFitMenu &CPFitOpt)
  {
    int a;
    is >> a;
    CPFitOpt = static_cast<CPFitMenu>(a);

    return is;
  }

  inline std::istream &operator>>(std::istream &is, GenVars &GenVarsOpt)
  {
    int a;
    is >> a;
    GenVarsOpt = static_cast<GenVars>(a);

    return is;
  }

  inline std::istream &operator>>(std::istream &is, OmegaRec &OmegaRecOpt)
  {
    int a;
    is >> a;
    OmegaRecOpt = static_cast<OmegaRec>(a);

    return is;
  }

  inline std::istream &operator>>(std::istream &is, Efficiency &EfficiencyOpt)
  {
    int a;
    is >> a;
    EfficiencyOpt = static_cast<Efficiency>(a);

    return is;
  }

  inline std::istream &operator>>(std::istream &is, DataType &DataTypeOpt)
  {
    int a;
    is >> a;
    DataTypeOpt = static_cast<DataType>(a);

    return is;
  }

  class Menu
  {
  private:
    std::map<int, TString> MenuOpt;
    std::vector<TString> fileOpt;

    const int ChooseMenu;
    const TString
        MenuName[8] = {"Main Menu", "NeutRec Menu", "Data Type", "Analysis file", "Final CPV Fit", "Generated Variables", "Omega-pi0 Reconstruction Menu", "Efficiency Check Menu"},

        ChooseOpt = "Choose the option: ";

  public:
    Menu(int ChooseMenu) : ChooseMenu(ChooseMenu)
    {
      switch (ChooseMenu)
      {
      case 0:
      {
        break;
      }
      case 1:
      {
        MenuOpt[int(NeutRecMenu::BARE_TRILATERATION)] = Form("%d. Bare trilateration.", int(NeutRecMenu::BARE_TRILATERATION));
        MenuOpt[int(NeutRecMenu::TRILATERATION_KIN_FIT)] = Form("%d. Trilateration with kinematic fit.", int(NeutRecMenu::TRILATERATION_KIN_FIT));
        MenuOpt[int(NeutRecMenu::TRIANGLE)] = Form("%d. Trilateration with triangle.", int(NeutRecMenu::TRIANGLE));
        MenuOpt[int(NeutRecMenu::COMP_OF_MET)] = Form("%d. Comparison of methods, plotting, etc.", int(NeutRecMenu::COMP_OF_MET));
        MenuOpt[int(NeutRecMenu::EXIT)] = Form("%d. Exit.", int(NeutRecMenu::EXIT));

        break;
      }
      case 2:
      {
        MenuOpt[int(DataType::SIGNAL_TOT)] = Form("%d. Total MC signal.", int(DataType::SIGNAL_TOT));
        MenuOpt[int(DataType::SIG_BCG)] = Form("%d. Signal + background MC.", int(DataType::SIG_BCG));
        MenuOpt[int(DataType::MC_DATA)] = Form("%d. MC + Data.", int(DataType::MC_DATA));

        break;
      }
      case 4:
      {
        MenuOpt[int(CPFitMenu::HALF_SIGNAL_MC)] = Form("%d. Only Signal; 1/2 MC vs. 1/2 'Data'", int(CPFitMenu::HALF_SIGNAL_MC));
        MenuOpt[int(CPFitMenu::HALF_SIG_BCG_MC)] = Form("%d. Signal + Bcg; 1/2 MC vs. 1/2 'Data'", int(CPFitMenu::HALF_SIG_BCG_MC));
        MenuOpt[int(CPFitMenu::MC_DATA)] = Form("%d. Signal + Bcg; MC vs. Data", int(CPFitMenu::MC_DATA));
        MenuOpt[int(CPFitMenu::EXIT)] = Form("%d. Exit.", int(CPFitMenu::EXIT));

        break;
      }
      case 5:
      {
        MenuOpt[int(GenVars::GEN_VARS)] = Form("%d. Genenerated variables to channel mapper.", int(GenVars::GEN_VARS));
        MenuOpt[int(GenVars::SPLIT_CHANN)] = Form("%d. Split stream into decay channels.", int(GenVars::SPLIT_CHANN));
        MenuOpt[int(CPFitMenu::EXIT)] = Form("%d. Exit.", int(CPFitMenu::EXIT));

        break;
      }
      case 6:
      {
        MenuOpt[int(OmegaRec::OMEGA_REC)] = Form("%d. Reconstruction of Omega channel based on the kinematic fit.", int(OmegaRec::OMEGA_REC));
        MenuOpt[int(OmegaRec::OMEGA_CUTS)] = Form("%d. Adjustment of cuts to reject Omega channel.", int(OmegaRec::OMEGA_CUTS));
        MenuOpt[int(OmegaRec::PLOTS)] = Form("%d. Plots, comparisons, etc.", int(OmegaRec::PLOTS));
        MenuOpt[int(OmegaRec::EXIT)] = Form("%d. Exit.", int(OmegaRec::EXIT));

        break;
      }
      case 7:
      {
        MenuOpt[int(Efficiency::EFF_SCAN)] = Form("%d. Scan of the efficiency for the range of chosen cuts.", int(Efficiency::EFF_SCAN));
        MenuOpt[int(Efficiency::EFF_DIST)] = Form("%d. Distribution of the efficiency for the range of chosen cuts.", int(Efficiency::EFF_DIST));

        MenuOpt[int(Efficiency::EXIT)] = Form("%d. Exit.", int(Efficiency::EXIT));

        break;
      }
      }
    };

    Menu(int ChooseMenu, std::vector<TString> fileNames) : ChooseMenu(ChooseMenu), fileOpt(fileNames)
    {
      switch (ChooseMenu)
      {
      case 0:
      {
        break;
      }
      case 1:
      {
        MenuOpt[int(NeutRecMenu::BARE_TRILATERATION)] = Form("%d. Bare trilateration.", int(NeutRecMenu::BARE_TRILATERATION));
        MenuOpt[int(NeutRecMenu::TRILATERATION_KIN_FIT)] = Form("%d. Trilateration with kinematic fit.", int(NeutRecMenu::TRILATERATION_KIN_FIT));
        MenuOpt[int(NeutRecMenu::TRIANGLE)] = Form("%d. Trilateration with triangle.", int(NeutRecMenu::TRIANGLE));
        MenuOpt[int(NeutRecMenu::COMP_OF_MET)] = Form("%d. Comparison of methods, plotting, etc.", int(NeutRecMenu::COMP_OF_MET));
        MenuOpt[int(NeutRecMenu::EXIT)] = Form("%d. Exit.", int(NeutRecMenu::EXIT));

        break;
      }
      case 2:
      {
        MenuOpt[int(DataType::SIGNAL_TOT)] = Form("%d. Total MC signal.", int(DataType::SIGNAL_TOT));
        MenuOpt[int(DataType::SIG_BCG)] = Form("%d. Signal + background MC.", int(DataType::SIG_BCG));
        MenuOpt[int(DataType::MC_DATA)] = Form("%d. MC + Data.", int(DataType::MC_DATA));

        break;
      }
      case 4:
      {
        MenuOpt[int(CPFitMenu::HALF_SIGNAL_MC)] = Form("%d. Only Signal; 1/2 MC vs. 1/2 'Data'", int(CPFitMenu::HALF_SIGNAL_MC));
        MenuOpt[int(CPFitMenu::HALF_SIG_BCG_MC)] = Form("%d. Signal + Bcg; 1/2 MC vs. 1/2 'Data'", int(CPFitMenu::HALF_SIG_BCG_MC));
        MenuOpt[int(CPFitMenu::MC_DATA)] = Form("%d. Signal + Bcg; MC vs. Data", int(CPFitMenu::MC_DATA));
        MenuOpt[int(CPFitMenu::EXIT)] = Form("%d. Exit.", int(CPFitMenu::EXIT));

        break;
      }
      case 5:
      {
        MenuOpt[int(GenVars::GEN_VARS)] = Form("%d. Genenerated variables to channel mapper.", int(GenVars::GEN_VARS));
        MenuOpt[int(GenVars::SPLIT_CHANN)] = Form("%d. Split stream into decay channels.", int(GenVars::SPLIT_CHANN));
        MenuOpt[int(GenVars::EXIT)] = Form("%d. Exit.", int(GenVars::EXIT));

        break;
      }
      case 6:
      {
        MenuOpt[int(OmegaRec::OMEGA_REC)] = Form("%d. Reconstruction of Omega channel based on the kinematic fit.", int(OmegaRec::OMEGA_REC));
        MenuOpt[int(OmegaRec::OMEGA_CUTS)] = Form("%d. Adjustment of cuts to reject Omega channel.", int(OmegaRec::OMEGA_CUTS));
        MenuOpt[int(OmegaRec::PLOTS)] = Form("%d. Plots, comparisons, etc.", int(OmegaRec::PLOTS));
        MenuOpt[int(OmegaRec::EXIT)] = Form("%d. Exit.", int(OmegaRec::EXIT));

        break;
      }
      case 7:
      {
        MenuOpt[int(Efficiency::EFF_SCAN)] = Form("%d. Scan of the efficiency for the range of chosen cuts.", int(Efficiency::EFF_SCAN));
        MenuOpt[int(Efficiency::EFF_DIST)] = Form("%d. Distribution of the efficiency for the range of chosen cuts.", int(Efficiency::EFF_DIST));

        MenuOpt[int(Efficiency::EXIT)] = Form("%d. Exit.", int(Efficiency::EXIT));

        break;
      }
      }
    };

    void InitMenu() { std::cout << MenuName[ChooseMenu] << std::endl; };
    void EndMenu() { std::cout << ChooseOpt; };
    void ShowOpt()
    {
      switch (ChooseMenu)
      {
      case 0:
      {
        break;
      }
      case 1:
      {
        for (int i = 1; i < int(NeutRecMenu::OPT_TOT); i++)
          std::cout << MenuOpt[i] << std::endl;

        break;
      }
      case 2:
      {
        for (int i = 1; i < int(DataType::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 3:
      {
        int iter = 0;
        for (const TString &fileName : fileOpt)
        {
          std::cout << iter << ". " << fileName << std::endl;
          iter++;
        }
      }
      case 4:
      {
        for (int i = 1; i < int(CPFitMenu::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 5:
      {
        for (int i = 1; i < int(GenVars::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 6:
      {
        for (int i = 1; i < int(OmegaRec::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 7:
      {
        for (int i = 1; i < int(Efficiency::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      }
    }
  };
}

#endif