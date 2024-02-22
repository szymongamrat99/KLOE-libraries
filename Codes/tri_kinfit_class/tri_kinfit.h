#ifndef TRI_KINFIT_H
#define TRI_KINFIT_H

#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TError.h>

const Float_t TRF = 2.715; // ns

namespace KLOE
{
  class TrilaterationRec
  {
  private:
    Int_t _N_free, _N_const, _M, _jmin, _jmax, _loopcount;

    TMatrixD _V, _D, _D_T, _V_final, _V_aux, _V_min, _Aux, _V_invert, _V_init;

    TVectorD _X, _C, _X_final, _L, _CORR, _X_init, _X_min, _C_min, _L_min, _C_aux, _L_aux, _X_init_min, _X_init_aux;

  public:
    TrilaterationRec(UInt_t N_free, UInt_t N_const, UInt_t M, UInt_t j, UInt_t loopcount);
    ~TrilaterationRec();
  };

  TrilaterationRec::TrilaterationRec(UInt_t N_free, UInt_t N_const, UInt_t M, UInt_t j, UInt_t loopcount)
  {
    _N_free = N_free;
    _N_const = N_const;
    _M = M;

    _jmin = -j;
    _jmax = j;

    _loopcount = loopcount;

    _V.ResizeTo(_N_free + _N_const, _N_free + _N_const);
    _D.ResizeTo(_M, _N_free + _N_const);
    _D_T.ResizeTo(_)
  }

  TrilaterationRec::~TrilaterationRec()
  {
  }

}

#endif
