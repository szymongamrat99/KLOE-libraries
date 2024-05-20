#ifndef KINFIT_H
#define KINFIT_H

#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TError.h>
#include <TF1.h>

const Float_t TRF = 2.715; // ns

namespace KLOE
{
  class KinFit
  {
  private:
    Int_t _N_free, _N_const, _M, _jmin, _jmax, _loopcount;
    Float_t _CHISQR, _FUNVAL, _CHISQRTMP, _FUNVALTMP, _CHISQRMIN, _FUNVALMIN;
    Double_t _det;

    TMatrixD _V, _D, _D_T, _V_final, _V_aux, _V_min, _Aux, _V_invert, _V_init;

    TVectorD _X, _C, _X_final, _L, _CORR, _X_init, _X_min, _C_min, _L_min, _C_aux, _L_aux, _X_init_min, _X_init_aux;

    Bool_t _fail;

  public:
    KinFit(UInt_t N_free, UInt_t N_const, UInt_t M, UInt_t j, UInt_t loopcount);
    KinFit(UInt_t N_free, UInt_t N_const, UInt_t M, UInt_t loopcount);
    ~KinFit();

    Int_t ParameterInitialization(Float_t *Params, Float_t *Errors);

    Int_t FitFunction();

  };

}

#endif
