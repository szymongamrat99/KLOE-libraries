#if !defined(CONSTRAINTS_TRI_H)
#define CONSTRAINTS_TRI_H

#include "../../Include/Codes/neutral_mom.h"
#include "../../Include/Codes/lorentz_transf.h"
#include "../../Include/Codes/plane_intersection.h"
#include "../../Include/Codes/closest_approach.h"

#include "../../Include/const.h"

Double_t ene_consv(Double_t *, Double_t *);
Double_t minv_consv(Double_t *, Double_t *);
Double_t x_consv(Double_t *, Double_t *);
Double_t y_consv(Double_t *, Double_t *);
Double_t z_consv(Double_t *, Double_t *);


#endif // CONSTRAINTS_TRI_H
