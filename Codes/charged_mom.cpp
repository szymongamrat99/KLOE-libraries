#include "charged_mom.h"
#include "TMath.h"
#include "../const.h"

void charged_mom(Float_t Curv, Float_t Phiv, Float_t Cotv, Float_t *mom_vec, UInt_t mode = 1)
{
    mom_vec[0] = cos(Phiv)*1000./abs( Curv );
    mom_vec[1] = sin(Phiv)*1000./abs( Curv );
    mom_vec[2] = Cotv*1000./abs( Curv );

    if(mode == 1) mom_vec[3] = sqrt(pow(mom_vec[0],2) + pow(mom_vec[1],2) + pow(mom_vec[2],2) + pow(m_pich,2));
    else if(mode == 2) mom_vec[3] = sqrt(pow(mom_vec[0],2) + pow(mom_vec[1],2) + pow(mom_vec[2],2) + pow(m_ele,2));
    else if(mode == 3) mom_vec[3] = sqrt(pow(mom_vec[0],2) + pow(mom_vec[1],2) + pow(mom_vec[2],2) + pow(m_mu,2));
}