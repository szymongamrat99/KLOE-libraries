#include "triple_gaus.h"

Double_t triple_gaus(Double_t *x, Double_t *p)
{
  Double_t norm1 = 1./(p[2]*sqrt(2*M_PI));
  Double_t norm2 = 1./(p[5]*sqrt(2*M_PI));
  Double_t norm3 = 1./(p[8]*sqrt(2*M_PI));

  Double_t gaus1 = p[0]*norm1*exp(-pow((x[0] - p[1])/p[2],2)/2.);
  Double_t gaus2 = p[3]*norm2*exp(-pow((x[0] - p[4])/p[5],2)/2.);
  Double_t gaus3 = p[6]*norm3*exp(-pow((x[0] - p[7])/p[8],2)/2.);

  Double_t value = gaus1 + gaus2 + gaus3;

  return value;
}

Double_t comb_std_dev(const Double_t *p)
{
  Double_t nominator = p[0] * p[2] + p[3] * p[5] + p[6] * p[8];
  Double_t denominator = p[0] + p[3] + p[6];

  Double_t std_dev = nominator / denominator;

  return std_dev;
}