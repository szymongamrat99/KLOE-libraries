#ifndef CONSTRAINTSTRI_H
#define CONSTRAINTSTRI_H

#include <TMath.h>
#include <TF1.h>

#include "../../../const.h"
#include "../../reconstructor.h"
#include "../../charged_mom.h"
#include "../../neutral_mom.h"
#include "../../lorentz_transf.h"
#include "../inc/KinFit.h"

/*List of variables (order is the following)*/
/*
  For 4 gamma decay, without charged decay:
  4 x 5 (Xcl, Ycl, Zcl, Tcl, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, sqrt(S))

  IP determined during fit (using direction of Kne momentum from clusters)
  Neutral kaon 4-vec determined during fit

--------------------------------------------------------------------------------

  For 4 gamma decay, with charged decay:
  2 x 3 (Curv, Phiv, Cotv)
  4 x 5 (Xcl, Ycl, Zcl, Tcl, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, sqrt(S))

  IP determined during fit (using direction of Kch momentum)
  Charged kaon 4-mom determined during fit (using the boost method as well)
  Neutral kaon 4-vec determined during fit
  Neutral kaon dependent of Kch direction and energy

--------------------------------------------------------------------------------

  For 6 gamma decay, without charged decay:
  6 x 5 (Xcl, Ycl, Zcl, Tcl, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, sqrt(S))

  IP determined during fit (using direction of Kne momentum from clusters)
  Neutral kaon 4-vec determined during fit

--------------------------------------------------------------------------------

  For 6 gamma decay, with charged decay:
  2 x 3 (Curv, Phiv, Cotv)
  6 x 5 (Xcl, Ycl, Zcl, Tcl, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, sqrt(S))

  IP determined during fit (using direction of Kch momentum)
  Charged kaon 4-mom determined during fit (using the boost method as well)
  Neutral kaon 4-vec determined during fit
  Neutral kaon dependent of Kch direction and energy

--------------------------------------------------------------------------------

*/

namespace KLOE
{
  class KinFitSignal : public KinFit
  {
  private:
    std::vector<Int_t> _chosen_constraints;

    std::vector<ChPart> _PiCh;
    std::vector<NeuPart> _Photon;
    std::vector<NeuPart> _PiNeu;
    std::vector<Kaon> _Kaon;
    Phi _PhiMeson;

    std::vector<Double_t> p;
    const std::vector<Double_t> p_const;

    using func = Double_t (KinFitSignal::*)(Double_t *, Double_t *);
    std::map<Int_t, func> _this_constraints;

    std::map<Int_t, TString> _constraint_name;

  public:
    KinFitSignal(Int_t N_free, Int_t N_const, Int_t N_clus, Int_t M, std::vector<Double_t> p_init);

    /* Specific physical constraints */

    // Pairing of photons to pi0
    void PhotonPairing();

    //
    // Invariant mass conservation - trilateration
    Double_t MinvConsvNeu(Double_t *, Double_t *);
    Double_t MinvConsvCh(Double_t *, Double_t *);

    // Energy conservation in LAB
    Double_t EneConsvLAB(Double_t *, Double_t *);
    Double_t PxConsvLAB(Double_t *, Double_t *);
    Double_t PyConsvLAB(Double_t *, Double_t *);
    Double_t PzConsvLAB(Double_t *, Double_t *);

    // Gamma TOF
    Double_t Gamma1TOF(Double_t *, Double_t *);
    Double_t Gamma2TOF(Double_t *, Double_t *);
    Double_t Gamma3TOF(Double_t *, Double_t *);
    Double_t Gamma4TOF(Double_t *, Double_t *);

    Double_t AngleChMom(Double_t *, Double_t *);
  };

}

#endif
