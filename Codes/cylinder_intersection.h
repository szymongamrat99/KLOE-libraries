#ifndef CYLINDER_INTERSECTION_H
#define CYLINDER_INTERSECTION_H

#include <iostream>
#include <TMath.h>

#define Rmax 200
#define zmax 165


Int_t inter_point(Float_t *, Float_t *, Float_t *);
void barrel_inter(Float_t *, Float_t *, Float_t *);
void endcap_inter(Float_t *, Float_t *, Float_t *);

#endif // !CYLINDER_INTERSECTION_H