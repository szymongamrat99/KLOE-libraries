#if !defined(CONSTRAINTS_OMEGA_H)
#define CONSTRAINTS_OMEGA_H

#include <neutral_mom.h>
#include <lorentz_transf.h>
#include <plane_intersection.h>
#include <closest_approach.h>
#include <pi0_photon_pair.h>

#include <reconstructor.h>

#include <const.h>

namespace OmegaConstraints
{

  Double_t ene_consv(Double_t *, Double_t *);
  Double_t px_consv(Double_t *, Double_t *);
  Double_t py_consv(Double_t *, Double_t *);
  Double_t pz_consv(Double_t *, Double_t *);
  Double_t minv_omega_consv(Double_t *, Double_t *);
  Double_t minv_pi01_consv(Double_t *, Double_t *);
  Double_t minv_pi02_consv(Double_t *, Double_t *);
  Double_t x_consv(Double_t *, Double_t *);
  Double_t y_consv(Double_t *, Double_t *);
  Double_t z_consv(Double_t *, Double_t *);

};

#endif // CONSTRAINTS_OMEGA_H
