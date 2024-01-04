#ifndef LINE_EMC_INTERSECTION_H
#define LINE_EMC_INTERSECTION_H
#define RT 200    //cm - EMC radius
#define zmax 165  //cm - max z length
#define zmin -165 //cm - min z length

#include <iostream>
#include <TMath.h>
#include <TVector3.h>

Int_t cluster_finder(Float_t *, Float_t *, Float_t *);

#endif