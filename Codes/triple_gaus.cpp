#include "triple_gaus.h"
#include <iostream>

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
  Double_t mean[3] = {p[1], p[4], p[7]}, sigma[3] = {p[2], p[5], p[8]}, N[3] = {p[0], p[3], p[6]};

  Double_t NTot = N[0] + N[1] + N[2];

  Double_t nominator = sqrt(pow(N[0]*sigma[0],2) + pow(N[1]*sigma[1],2) + pow(N[2]*sigma[2],2));

  Double_t std_dev = nominator / NTot;

  return std_dev;
}

Double_t comb_std_dev_err(const Double_t *p, const Double_t *err)
{
  Double_t mean[3] = {p[1], p[4], p[7]}, sigma[3] = {p[2], p[5], p[8]}, N[3] = {p[0], p[3], p[6]};

  Double_t mean_err[3] = {err[1], err[4], err[7]}, sigma_err[3] = {err[2], err[5], err[8]}, N_err[3] = {err[0], err[3], err[6]};

  Double_t NTot = N[0] + N[1] + N[2];

  Double_t nominator = sqrt(pow(N[0]*sigma[0],2) + pow(N[1]*sigma[1],2) + pow(N[2]*sigma[2],2));

  Double_t err_denominator = pow(NTot,2)*nominator;
  Double_t err_sigma_tot = 0., err_N_tot = 0., err_nominator = 0.;

  for(Int_t i = 0; i < 3; i++)
  {
    err_sigma_tot += pow(sigma_err[i],2)*pow(N[i],4)*pow(NTot,2)*pow(sigma[i],2);
  }

  err_N_tot += pow(N_err[0],2)*pow(pow(N[1]*sigma[1],2)+pow(N[2]*sigma[2],2)-N[0]*(N[1]+N[2])*pow(sigma[0],2),2);

  err_N_tot += pow(N_err[1],2)*pow(pow(N[0]*sigma[0],2)+pow(N[2]*sigma[2],2)-N[1]*(N[0]+N[2])*pow(sigma[1],2),2);

  err_N_tot += pow(N_err[2],2)*pow(pow(N[0]*sigma[0],2)+pow(N[1]*sigma[1],2)-N[2]*(N[0]+N[1])*pow(sigma[2],2),2);

  err_nominator = err_sigma_tot + err_N_tot;


  Double_t std_dev_err = sqrt(err_nominator) / err_denominator;

  return std_dev_err;
}