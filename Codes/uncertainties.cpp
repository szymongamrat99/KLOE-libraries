#include "uncertainties.h"
#include "TMath.h"
#include "../const.h"

Double_t clu_ene_error(Float_t Enecl)
{
    Double_t sigma_E = 0.;

    //Error in MeV

    sigma_E = Enecl*0.057/sqrt( 0.001*Enecl );

    if( Enecl > 0. ) return sigma_E;
    else return -999;
};

Double_t clu_time_error(Float_t Enecl)
{
    Double_t sigma_time = 0.;

    //Error in ns

    sigma_time = sqrt(pow(0.054/sqrt( 0.001*Enecl ),2) + pow(0.05,2) );

    if( Enecl > 0. ) return sigma_time;
    else return -999;
}