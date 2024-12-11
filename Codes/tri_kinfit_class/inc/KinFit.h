#ifndef KINFIT_H
#define KINFIT_H

#include <vector>
#include<algorithm>
#include <utility> // dla std::swap

#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TError.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include "../../../const.h"

namespace KLOE
{
  class KinFit
  {
  private:
    Float_t _CHISQR, _FUNVAL, _CHISQRTMP, _FUNVALTMP, _CHISQRMIN, _FUNVALMIN;
    Double_t _det;

    TMatrixD _V, _D, _D_T, _V_final, _V_aux, _V_min, _Aux, _V_invert, _V_init;

    TVectorD _X, _C, _X_final, _L, _CORR, _X_init, _X_min, _C_min, _L_min, _C_aux, _L_aux, _X_init_min, _X_init_aux; 

  protected:
    const Double_t _TRF = 2.715; // ns
    Int_t _N_free, _N_const, _M, _M_act;
    Bool_t _err_flag, _fail;
    Int_t *_selected, _counter = 0, _jmin, _jmax, _loopcount, _N_clus;
    std::vector<TF1*> _constraints;
    Double_t _value_min;

    std::vector<Int_t> _chosen;

    struct NeuPart
    {
      const UInt_t index;
      Int_t pid; // Geant4 pids
      TLorentzVector FourMom, FourPos;
      TVector3 boost;

      Double_t cluster[5];
      Double_t InvMass;
    };

    struct ChPart
    {
      Int_t pid; // Geant4 pids
      TLorentzVector FourMom, FourPos;
      TVector3 boost;

      Double_t TrkParameters[3];
      Double_t InvMass;
    };

    struct Kaon
    {
      Int_t pid; // Geant4 pids
      TLorentzVector FourMom, FourPos;
      TVector3 boost;

      Double_t InvMass;
    };

    struct Phi
    {
      Int_t pid; // Geant4 pids
      TLorentzVector FourMom, FourPos;
      TVector3 boost;

      Double_t InvMass;
    };
    



  public:
    KinFit();
    ~KinFit();

    Int_t ParameterInitialization(Float_t *Params, Float_t *Errors);

    Int_t ConstraintSet(std::vector<TF1*> ConstrSet);

    void PhotonPairing(std::vector<NeuPart> _Photons);

    void FitFunction();

    Double_t EnergyCalc(Double_t *p, Double_t mass);
    Double_t EnergyCalc(TLorentzVector p, Double_t mass);
  };

}

#endif
