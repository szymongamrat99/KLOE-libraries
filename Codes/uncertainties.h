#ifndef UNCERTAINTIES_H
#define UNCERTAINTIES_H

#include "TMath.h"
#include "../const.h"

Double_t clu_ene_error(Float_t);
Double_t clu_time_error(Float_t);
Double_t clu_x_error();
Double_t clu_y_error();
Double_t clu_z_error(Float_t);

#endif