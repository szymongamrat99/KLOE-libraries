#ifndef MAIN_MENU_H
#define MAIN_MENU_H

#include <string>
#include <map>

namespace Controls
{
  enum class MainMenu
  {
    GEN_VARS = 1,
    KCHREC_NO_BOOST = 2,
    KCHREC_BOOST = 3,
    IP_EV_BY_EV = 4,
    KNEREC_TRILAT = 5,
    KNEREC_TRIANGLE = 6,
    EFF_SIG_TO_BCG = 7,
    KIN_FITS = 8,
    TRANSF_TO_CM = 9,
    CPV_NORM = 10,
    EXIT = 11
  };

  std::istream &operator>>(std::istream &is, MainMenu &mainMenuOpt)
  {
    int a;
    is >> a;
    mainMenuOpt = static_cast<MainMenu>(a);

    return is;
  }

  class Menu
  {
  private:
    std::map<MainMenu, TString> MenuOpt;
    const int ChooseMenu, MenuNum = 5;
    const TString
        MenuName[5] = {"Main Menu", "", "", "", ""},
        
        ChooseOpt = "Choose the option: ";

  public:
    Menu(int ChooseMenu) : ChooseMenu(ChooseMenu) 
    {
                            
    };

    void InitMenu() { std::cout << MenuName[ChooseMenu] << std::endl; };
    void EndMenu() { std::cout << ChooseOpt; };
  };
}

#endif