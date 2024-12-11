// File with methods from KinFit class for KLOE
// Date: 05/03/2024
// Author: Szymon Gamrat

#include "../inc/KinFit.h"

using namespace KLOE;

Int_t KinFit::ParameterInitialization(Float_t *Params, Float_t *Errors)
{
  // Passed arrays should be of size _N_free + _N_const

  for (Int_t i = 0; i < _N_free + _N_const; i++)
  {
    _X(i) = Params[i];

    if (i < _N_free)
    {
      _V(i, i) = pow(Errors[i], 2);
      if (std::isnan(_V(i, i)))
        _err_flag = true;
    }
  }

  return 0;
};

void KinFit::FitFunction()
{
  for (Int_t l = 0; l < _M; l++)
  {
    _C(l) = _constraints[l]->EvalPar(0, _X.GetMatrixArray());
    for (Int_t m = 0; m < _N_free + _N_const; m++)
    {
      _constraints[l]->SetParameters(_X.GetMatrixArray());
      if (m < _N_free)
        _D(l, m) = _constraints[l]->GradientPar(m, 0, 0.01);
      else
        _D(l, m) = 0;
    }
  }

  _D_T.Transpose(_D);

  _Aux = (_D * _V * _D_T);

  _Aux.Invert(&_det);

  if (!TMath::IsNaN(_det) && _det != 0)
  {
    _L = (_Aux * _C);

    _CORR = _V * _D_T * _L;

    _X_final = _X - _CORR;
    _V_final = _V - _V * _D_T * _Aux * _D * _V;

    for (Int_t l = 0; l < _M; l++)
      _C(l) = _constraints[l]->EvalPar(0, _X_final.GetMatrixArray());

    _FUNVAL = Dot((_X_final - _X_init).GetSub(0, _N_free - 1), _V_invert * (_X_final - _X_init).GetSub(0, _N_free - 1)) + Dot(_L, _C);

    _CHISQR = Dot((_X_final - _X_init).GetSub(0, _N_free - 1), _V_invert * (_X_final - _X_init).GetSub(0, _N_free - 1));
  }
  else
  {
    _fail = 1;
  }

  if (_fail == 0)
  {
    _X = _X_final;
    _X_init_aux = _X_init;
    _V_aux = _V_final;
    _L_aux = _L;
    _C_aux = _C;
    _FUNVALTMP = _FUNVAL;
    _CHISQRTMP = _CHISQR;
  }
};

Double_t KinFit::EnergyCalc(Double_t *p, Double_t mass)
{
  Double_t energy = 0.;

  for (Int_t i = 0; i < 3; i++)
    energy += pow(p[i], 2);

  energy += pow(mass, 2);
  energy = sqrt(energy);

  return energy;
};

Double_t KinFit::EnergyCalc(TLorentzVector p, Double_t mass)
{
  Double_t energy = 0.;

  for (Int_t i = 0; i < 3; i++)
    energy += pow(p[i], 2);

  energy += pow(mass, 2);
  energy = sqrt(energy);

  return energy;
};

void KinFit::PhotonPairing(std::vector<NeuPart> _Photons)
{
  Int_t PhotonsNum = _Photons.size();

  if (PhotonsNum % 2 != 0)
  {
    std::cerr << "Liczba fotonów musi być parzysta." << std::endl;
    return;
  }

  const Int_t div = Int_t(PhotonsNum / 2.);
  Double_t AuxMinv[div];
  Double_t AuxChi2 = 0., AuxChi2Min = 99999.;

  do
  {
    for (Int_t i = 0; i < div; i++)
    {
      AuxMinv[i] += pow(_Photons[i * div].FourMom[3] + _Photons[i * div + 1].FourMom[3], 2);

      for (Int_t j = 0; j < 3; j++)
      {
        AuxMinv[i] -= pow(_Photons[i * div].FourMom[j] + _Photons[i * div + 1].FourMom[j], 2);
      }

      AuxMinv[i] = sqrt(AuxMinv[i]);

      AuxChi2 += pow(AuxMinv[i] - mPi0, 2);
    }

    AuxChi2 = sqrt(AuxChi2);

    if (AuxChi2 < AuxChi2Min)
    {
      AuxChi2Min = AuxChi2;
    }

  } while (std::next_permutation(_Photons.begin(), _Photons.end(), [](const auto & lhs, const auto & rhs) 
                                 { return lhs.index < rhs.index; }));
}