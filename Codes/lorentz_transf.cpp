#include "lorentz_transf.h"

void lorentz_transf(Float_t *boost_vec, Float_t *vec_init, Float_t *vec_end)
{
  TVector3 boost(boost_vec);
  TLorentzVector before(vec_init);

  before.Boost(boost);

  for(Int_t i = 0; i < 4; i++)
    vec_end[i] = before(i);
}