#include <string.h>

#include "constraints_tri.h"

const Double_t Trf = 2.715;

Double_t ene_consv(Double_t *x, Double_t *p)
{
  Float_t boost_vec[3] = {-p[20]/p[23], -p[21]/p[23], -p[22]/p[23]};
  Float_t gamma_mom[4][4], vec_init[4], vec_end[4], neu_vtx[4] = {p[24], p[25], p[26], p[27]}; 

  neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
  neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
  neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
  neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

  vec_init[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
  vec_init[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
  vec_init[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
  vec_init[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

  lorentz_transf(boost_vec, vec_init, vec_end);

  Double_t value = vec_end[3] - (p[23]/2.);
  return value;
}

Double_t minv_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], neu_vtx[4] = {p[24], p[25], p[26], p[27]}, inv_mass_kaon; 

  neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
  neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
  neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
  neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

  kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
  kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
  kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
  kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

  inv_mass_kaon = sqrt(pow(kaon_mom[3],2) - pow(kaon_mom[0],2) - pow(kaon_mom[1],2) - pow(kaon_mom[2],2));

  Double_t value = inv_mass_kaon - m_k0;
  return value;
}

Double_t x_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], neu_vtx[4] = {p[24], p[25], p[26], p[27]}, inv_mass_kaon,
          bhabha_vtx[3] = {p[28], p[29], p[30]}, y_axis[3] = {0., 1., 0.}, ip[3]; 


  neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
  neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
  neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
  neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

  kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
  kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
  kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
  kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

  plane_intersection(bhabha_vtx, y_axis, neu_vtx, kaon_mom, ip);

  ip[0] = p[28];
  ip[1] = p[29];
  if( abs(p[30] - ip[2]) > 2 ) ip[2] = p[30];

  Double_t kaon_vel[3] = {c_vel * kaon_mom[0]/kaon_mom[3], c_vel * kaon_mom[1]/kaon_mom[3], c_vel * kaon_mom[2]/kaon_mom[3]}, dist[3] = {neu_vtx[0] - ip[0], neu_vtx[1] - ip[1], neu_vtx[2] - ip[2]};

  Double_t kaon_vel_tot = sqrt(pow(kaon_vel[0],2) + pow(kaon_vel[1],2) + pow(kaon_vel[2],2)), kaon_path = sqrt(pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2));

  Double_t time_diff = neu_vtx[3] - (kaon_path/kaon_vel_tot);

  if(time_diff < -Trf/2. && time_diff > -3.*Trf/2.) neu_vtx[3] += Trf;
  if(time_diff < -3.*Trf/2. && time_diff > -5.*Trf/2.) neu_vtx[3] += 2*Trf;
  if(time_diff < -5.*Trf/2. && time_diff > -7.*Trf/2.) neu_vtx[3] += 3*Trf;
  if(time_diff < -9.*Trf/2. && time_diff > -9.*Trf/2.) neu_vtx[3] += 4*Trf;
  if(time_diff < -9.*Trf/2. && time_diff > -11.*Trf/2.) neu_vtx[3] += 5*Trf;
  if(time_diff > Trf/2. && time_diff < 3.*Trf/2.) neu_vtx[3] -= Trf;
  if(time_diff > 3.*Trf/2. && time_diff < 5.*Trf/2.) neu_vtx[3] -= 2*Trf;
  if(time_diff > 5.*Trf/2. && time_diff < 7.*Trf/2.) neu_vtx[3] -= 3*Trf;

  Double_t value = kaon_vel[0] * neu_vtx[3] - dist[0];

  return value;
}

Double_t y_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], neu_vtx[4] = {p[24], p[25], p[26], p[27]}, inv_mass_kaon,
          bhabha_vtx[3] = {p[28], p[29], p[30]}, y_axis[3] = {0., p[21], 0.}, ip[3]; 

  neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
  neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
  neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
  neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

  kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
  kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
  kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
  kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

  plane_intersection(bhabha_vtx, y_axis, neu_vtx, kaon_mom, ip);

  ip[0] = p[28];
  ip[1] = p[29];

  if( abs(p[30] - ip[2]) > 2 ) ip[2] = p[30];

  Double_t kaon_vel[3] = {c_vel * kaon_mom[0]/kaon_mom[3], c_vel * kaon_mom[1]/kaon_mom[3], c_vel * kaon_mom[2]/kaon_mom[3]}, dist[3] = {neu_vtx[0] - ip[0], neu_vtx[1] - ip[1], neu_vtx[2] - ip[2]};

  Double_t kaon_vel_tot = sqrt(pow(kaon_vel[0],2) + pow(kaon_vel[1],2) + pow(kaon_vel[2],2)), kaon_path = sqrt(pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2));

  Double_t time_diff = neu_vtx[3] - (kaon_path/kaon_vel_tot);

  if(time_diff < -Trf/2. && time_diff > -3.*Trf/2.) neu_vtx[3] += Trf;
  if(time_diff < -3.*Trf/2. && time_diff > -5.*Trf/2.) neu_vtx[3] += 2*Trf;
  if(time_diff < -5.*Trf/2. && time_diff > -7.*Trf/2.) neu_vtx[3] += 3*Trf;
  if(time_diff < -9.*Trf/2. && time_diff > -9.*Trf/2.) neu_vtx[3] += 4*Trf;
  if(time_diff < -9.*Trf/2. && time_diff > -11.*Trf/2.) neu_vtx[3] += 5*Trf;
  if(time_diff > Trf/2. && time_diff < 3.*Trf/2.) neu_vtx[3] -= Trf;
  if(time_diff > 3.*Trf/2. && time_diff < 5.*Trf/2.) neu_vtx[3] -= 2*Trf;
  if(time_diff > 5.*Trf/2. && time_diff < 7.*Trf/2.) neu_vtx[3] -= 3*Trf;

  Double_t value = kaon_vel[1] * neu_vtx[3] - dist[1];
  return value;
}

Double_t z_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], neu_vtx[4] = {p[24], p[25], p[26], p[27]}, inv_mass_kaon,
          bhabha_vtx[3] = {p[28], p[29], p[30]}, y_axis[3] = {0., p[21], 0.}, ip[3]; 

  neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
  neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
  neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
  neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

  kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
  kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
  kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
  kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

  plane_intersection(bhabha_vtx, y_axis, neu_vtx, kaon_mom, ip);

  ip[0] = p[28];
  ip[1] = p[29];

  if( abs(p[30] - ip[2]) > 2 ) ip[2] = p[30];

  Double_t kaon_vel[3] = {c_vel * kaon_mom[0]/kaon_mom[3], c_vel * kaon_mom[1]/kaon_mom[3], c_vel * kaon_mom[2]/kaon_mom[3]}, dist[3] = {neu_vtx[0] - ip[0], neu_vtx[1] - ip[1], neu_vtx[2] - ip[2]};

  Double_t kaon_vel_tot = sqrt(pow(kaon_vel[0],2) + pow(kaon_vel[1],2) + pow(kaon_vel[2],2)), kaon_path = sqrt(pow(dist[0],2) + pow(dist[1],2) + pow(dist[2],2));

  Double_t time_diff = neu_vtx[3] - (kaon_path/kaon_vel_tot);

  if(time_diff < -Trf/2. && time_diff > -3.*Trf/2.) neu_vtx[3] += Trf;
  if(time_diff < -3.*Trf/2. && time_diff > -5.*Trf/2.) neu_vtx[3] += 2*Trf;
  if(time_diff < -5.*Trf/2. && time_diff > -7.*Trf/2.) neu_vtx[3] += 3*Trf;
  if(time_diff < -9.*Trf/2. && time_diff > -9.*Trf/2.) neu_vtx[3] += 4*Trf;
  if(time_diff < -9.*Trf/2. && time_diff > -11.*Trf/2.) neu_vtx[3] += 5*Trf;
  if(time_diff > Trf/2. && time_diff < 3.*Trf/2.) neu_vtx[3] -= Trf;
  if(time_diff > 3.*Trf/2. && time_diff < 5.*Trf/2.) neu_vtx[3] -= 2*Trf;
  if(time_diff > 5.*Trf/2. && time_diff < 7.*Trf/2.) neu_vtx[3] -= 3*Trf;

  Double_t value = kaon_vel[2] * neu_vtx[3] - dist[2];
  return value;
}