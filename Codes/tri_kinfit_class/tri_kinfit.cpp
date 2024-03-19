// File with methods from TrilaterationRec class for KLOE
// Date: 05/03/2024
// Author: Szymon Gamrat

#include "tri_kinfit.h"

using namespace KLOE;

TrilaterationRec::TrilaterationRec(UInt_t N_free, UInt_t N_const, UInt_t M, UInt_t j, UInt_t loopcount)
{
  _N_free = N_free;
  _N_const = N_const;
  _M = M;

  _jmin = -j;
  _jmax = j;

  _loopcount = loopcount;

  // Definitions of all matrices and vectors

  _V.ResizeTo(_N_free + _N_const, _N_free + _N_const);

  _D.ResizeTo(_M, _N_free + _N_const);
  _D_T.ResizeTo(_N_free + _N_const, _M);

  _X.ResizeTo(_N_free + _N_const);
  _C.ResizeTo(_M);
  _L.ResizeTo(_M);
  _CORR.ResizeTo(_N_free + _N_const);
}
TrilaterationRec::TrilaterationRec(UInt_t N_free, UInt_t N_const, UInt_t M, UInt_t loopcount)
{
  _N_free = N_free;
  _N_const = N_const;
  _M = M;

  _loopcount = loopcount;

  // Definitions of all matrices and vectors

  _V.ResizeTo(_N_free + _N_const, _N_free + _N_const);

  _D.ResizeTo(_M, _N_free + _N_const);
  _D_T.ResizeTo(_N_free + _N_const, _M);

  _X.ResizeTo(_N_free + _N_const);
  _C.ResizeTo(_M);
  _L.ResizeTo(_M);
  _CORR.ResizeTo(_N_free + _N_const);
}
TrilaterationRec::~TrilaterationRec()
{
}

Int_t TrilaterationRec::ParameterInitialization(Float_t *Params, Float_t *Errors)
{
  // Passed arrays should be of size _N_free + _N_const

  for (Int_t i = 0; i < _N_free + _N_const; i++)
  {
    _X(i) = Params[i];

    if (i < _N_free)
    {
      _V(i, i) = pow(Errors[i], 2);
      if (std::isnan(_V(i, i)))
        err_flag = true;
    }
  }

  return 0;
};

Int_t TrilaterationRec::FitFunction()
{
  for (Int_t l = 0; l < _M; l++)
  {
    _C(l) = constraints[l]->EvalPar(0, _X.GetMatrixArray());
    for (Int_t m = 0; m < _N_free + _N_const; m++)
    {
      constraints[l]->SetParameters(_X.GetMatrixArray());
      if (m < _N_free)
        _D(l, m) = constraints[l]->GradientPar(m, 0, 0.01);
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
      _C(l) = constraints[l]->EvalPar(0, _X_final.GetMatrixArray());

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