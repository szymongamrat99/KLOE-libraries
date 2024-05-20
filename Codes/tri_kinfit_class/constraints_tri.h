#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TError.h>
#include <TF1.h>

const Float_t TRF = 2.715; // ns

namespace KLOE
{
  class Constraints
  {
  private:
    

  public:
    Constraints(UInt_t N_free, UInt_t N_const, UInt_t M, UInt_t j, UInt_t loopcount);
    Constraints(UInt_t N_free, UInt_t N_const, UInt_t M, UInt_t loopcount);
    ~Constraints();

    Double_t ParameterInitialization(Float_t *Params, Float_t *Errors);

    Int_t FitFunction();

  };

}

#endif
