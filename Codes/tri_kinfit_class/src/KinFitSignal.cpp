#include "../inc/KinFitSignal.h"
#include "../../../const.h"

using namespace KLOE;

// Standard constructor of the class

KinFitSignal::KinFitSignal(Int_t N_free, Int_t N_const, Int_t N_clus, Int_t M, std::vector<Double_t> p_init) : p(p_init), p_const(p_init)
{
  _N_free = N_free;
  _N_const = N_const;
  _N_clus = N_clus;
  _M = M;

  Int_t choice = 0;

  std::cout << "Which constraints to use?" << std::endl;
  std::cout << "1. Energy conservation" << std::endl;
  std::cout << "2. Px conservation" << std::endl;
  std::cout << "3. Py conservation" << std::endl;
  std::cout << "4. Pz conservation" << std::endl;
  std::cout << "5. Kch invariant mass" << std::endl;
  std::cout << "6. Kne invariant mass" << std::endl;
  std::cout << "7. 1st photon TOF triangle" << std::endl;
  std::cout << "8. 2nd photon TOF triangle" << std::endl;
  std::cout << "9. 3rd photon TOF triangle" << std::endl;
  std::cout << "10. 4rd photon TOF triangle" << std::endl;
  std::cout << "99. Enough" << std::endl;

  _this_constraints[1] = &KinFitSignal::EneConsvLAB;
  _this_constraints[2] = &KinFitSignal::PxConsvLAB;
  _this_constraints[3] = &KinFitSignal::PyConsvLAB;
  _this_constraints[4] = &KinFitSignal::PzConsvLAB;
  _this_constraints[5] = &KinFitSignal::MinvConsvCh;
  _this_constraints[6] = &KinFitSignal::MinvConsvNeu;
  _this_constraints[7] = &KinFitSignal::Gamma1TOF;
  _this_constraints[8] = &KinFitSignal::Gamma2TOF;
  _this_constraints[9] = &KinFitSignal::Gamma3TOF;
  _this_constraints[10] = &KinFitSignal::Gamma4TOF;

  _constraint_name[1] = "Energy_conservation";
  _constraint_name[2] = "Px_conservation";
  _constraint_name[3] = "Py_conservation";
  _constraint_name[4] = "Pz_conservation";
  _constraint_name[5] = "Minv_charged";
  _constraint_name[6] = "Minv_neutral";
  _constraint_name[7] = "Gamma1TOF";
  _constraint_name[8] = "Gamma2TOF";
  _constraint_name[9] = "Gamma3TOF";
  _constraint_name[10] = "Gamma4TOF";

  while(1)
  {
    if(_chosen_constraints.size() < _M)
    {
      std::cin >> choice;
      std::cout << "Choice: " << choice << std::endl;
    }
    else
      break;

    if(choice == 99)
      break;
    else
    {
      _constraints.push_back(new TF1(_constraint_name[choice],_this_constraints[choice], 0., 1., _N_free + _N_const));
    }
  }

  _M_act = _constraints.size();

  std::cout << "Total number of free parameters: " << _N_free << std::endl;
  std::cout << "Total number of const parameters: " << _N_const << std::endl;
  std::cout << "Total number of clusters: " << _N_clus << std::endl;
  std::cout << "Maximal number of constraints: " << _M << std::endl;
  std::cout << "Actual number of constraints: " << _M_act << std::endl;

};

Double_t KinFitSignal::MinvConsvCh(Double_t *x, Double_t *p)
{
  for (Int_t i = 0; i < 2; i++)
  {
    _PiCh[i].FourMom.SetPx(p[i * 2]);
    _PiCh[i].FourMom.SetPy(p[i * 2 + 1]);
    _PiCh[i].FourMom.SetPz(p[i * 2 + 2]);
    _PiCh[i].FourMom.SetE(EnergyCalc(_PiCh[i].FourMom, m_pich));
  }

  _Kaon[0].FourMom.SetPxPyPzE(_PiCh[0].FourMom[0] + _PiCh[1].FourMom[0],
                              _PiCh[0].FourMom[1] + _PiCh[1].FourMom[1],
                              _PiCh[0].FourMom[2] + _PiCh[1].FourMom[2],
                              _PiCh[0].FourMom[3] + _PiCh[1].FourMom[3]);

  _value_min = (_Kaon[0].FourMom.Mag() - m_k0);

  return _value_min;
};

/////////////////////////////////////////////////////////////////////////

Double_t KinFitSignal::AngleChMom(Double_t *x, Double_t *p)
{
  Double_t Kch_time;

  for (Int_t i = 0; i < 2; i++)
  {
    _PiCh[i].FourMom.SetPx(p[i * 2]);
    _PiCh[i].FourMom.SetPy(p[i * 2 + 1]);
    _PiCh[i].FourMom.SetPz(p[i * 2 + 2]);
    _PiCh[i].FourMom.SetE(EnergyCalc(_PiCh[i].FourMom, m_pich));
  }
  _Kaon[0].FourMom.SetPxPyPzE(_PiCh[0].FourMom[0] + _PiCh[1].FourMom[0],
                          _PiCh[0].FourMom[1] + _PiCh[1].FourMom[1],
                          _PiCh[0].FourMom[2] + _PiCh[1].FourMom[2],
                          _PiCh[0].FourMom[3] + _PiCh[1].FourMom[3]);

  _Kaon[0].FourPos.SetXYZT(p[6] - p[24], p[7] - p[25], p[8] - p[26], p[9] - p[27]);

  _PhiMeson.FourMom.SetPxPyPzE(p[20], p[21], p[22], p[23]);
  _PhiMeson.FourPos.SetXYZT(p[24], p[25], p[26], p[23]);

  _value_min = (_Kaon[0].FourMom.Mag() - m_k0);

  return _value_min;
};
