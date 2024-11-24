#include <constraints_omega.h>

namespace OmegaConstraints
{
  const Double_t Trf = 2.715; // ns
  Bool_t cond_detector = kFALSE;

  Int_t selected[4] = {1, 2, 3, 4};

  Double_t ene_consv(Double_t *x, Double_t *p)
  {
    Float_t
        gamma_mom[4][8],
        vec_init[4],
        vec_end[4],
        neu_vtx[4],
        trk[2][4];

    Double_t
        value[2],
        value_min,
        T0;

    T0 = p[30] * Trf;

    Reconstructor R;
    Solution S;

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3] + T0,
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }

    for (Int_t i = 0; i < 2; i++)
    {
      for (Int_t j = 0; j < 3; j++)
        trk[i][j] = p[20 + j + 3 * i];

      trk[i][3] = sqrt(pow(mPiCh, 2) + pow(trk[i][0], 2) + pow(trk[i][1], 2) + pow(trk[i][2], 2));
    }

    S = R.MySolve(selected);

    for (Int_t i = 0; i < 2; i++)
    {
      if (S.error[i] == 0)
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
        neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
        neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
        neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

        value[i] = trk[0][3] + trk[1][3] + gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3] - p[29];
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (abs(value[0]) < abs(value[1]))
      value_min = value[0];
    else
      value_min = value[1];

    return value_min;
  }

  Double_t px_consv(Double_t *x, Double_t *p)
  {
    Float_t
        gamma_mom[4][8],
        vec_init[4],
        vec_end[4],
        neu_vtx[4],
        trk[2][4];

    Double_t
        value[2],
        value_min,
        T0;

    T0 = p[30] * Trf;

    Reconstructor R;
    Solution S;

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3] + T0,
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }

    for (Int_t i = 0; i < 2; i++)
    {
      for (Int_t j = 0; j < 3; j++)
        trk[i][j] = p[20 + j + 3 * i];

      trk[i][3] = sqrt(pow(mPiCh, 2) + pow(trk[i][0], 2) + pow(trk[i][1], 2) + pow(trk[i][2], 2));
    }

    S = R.MySolve(selected);

    for (Int_t i = 0; i < 2; i++)
    {
      if (S.error[i] == 0)
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
        neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
        neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
        neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

        value[i] = trk[0][0] + trk[1][0] + gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0] - p[26];
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (abs(value[0]) < abs(value[1]))
      value_min = value[0];
    else
      value_min = value[1];

    return value_min;
  }

  Double_t py_consv(Double_t *x, Double_t *p)
  {
    Float_t
        gamma_mom[4][8],
        vec_init[4],
        vec_end[4],
        neu_vtx[4],
        trk[2][4];

    Double_t
        value[2],
        value_min,
        T0;

    T0 = p[30] * Trf;

    Reconstructor R;
    Solution S;

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3] + T0,
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }

    for (Int_t i = 0; i < 2; i++)
    {
      for (Int_t j = 0; j < 3; j++)
        trk[i][j] = p[20 + j + 3 * i];

      trk[i][3] = sqrt(pow(mPiCh, 2) + pow(trk[i][0], 2) + pow(trk[i][1], 2) + pow(trk[i][2], 2));
    }

    S = R.MySolve(selected);

    for (Int_t i = 0; i < 2; i++)
    {
      if (S.error[i] == 0)
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
        neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
        neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
        neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

        value[i] = trk[0][1] + trk[1][1] + gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1] - p[27];
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (abs(value[0]) < abs(value[1]))
      value_min = value[0];
    else
      value_min = value[1];

    return value_min;
  }

  Double_t pz_consv(Double_t *x, Double_t *p)
  {
    Float_t
        gamma_mom[4][8],
        vec_init[4],
        vec_end[4],
        neu_vtx[4],
        trk[2][4];

    Double_t
        value[2],
        value_min,
        T0;

    T0 = p[30] * Trf;

    Reconstructor R;
    Solution S;

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3] + T0,
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }

    for (Int_t i = 0; i < 2; i++)
    {
      for (Int_t j = 0; j < 3; j++)
        trk[i][j] = p[20 + j + 3 * i];

      trk[i][3] = sqrt(pow(mPiCh, 2) + pow(trk[i][0], 2) + pow(trk[i][1], 2) + pow(trk[i][2], 2));
    }

    S = R.MySolve(selected);

    for (Int_t i = 0; i < 2; i++)
    {
      if (S.error[i] == 0)
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
        neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
        neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
        neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

        value[i] = trk[0][2] + trk[1][2] + gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2] - p[28];
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (abs(value[0]) < abs(value[1]))
      value_min = value[0];
    else
      value_min = value[1];

    return value_min;
  }

  Double_t minv_omega_consv(Double_t *x, Double_t *p)
  {
    Int_t
        clusterind[4],
        clusterindpi0[2][2];

    Float_t
        gamma_mom[4][8],
        gamma_mom_pi0[2][2][4],
        pi0_mom[2][4],
        omega_mom[4],
        neu_vtx[4],
        trk[2][4],
        inv_mass_omega;

    Double_t value[2], value_min, T0;

    T0 = p[30] * Trf;

    Reconstructor R;
    Solution S;

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3] + T0,
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }

    S = R.MySolve(selected);

    for (Int_t i = 0; i < 2; i++)
    {
      for (Int_t j = 0; j < 3; j++)
        trk[i][j] = p[20 + j + 3 * i];

      trk[i][3] = sqrt(pow(mPiCh, 2) + pow(trk[i][0], 2) + pow(trk[i][1], 2) + pow(trk[i][2], 2));
    }

    for (Int_t i = 0; i < 2; i++)
    {

      if (!S.error[i])
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
        neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
        neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
        neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

        Pi0PhotonPair(clusterind, gamma_mom, clusterindpi0, gamma_mom_pi0, pi0_mom, true, trk, omega_mom);

        inv_mass_omega = sqrt(pow(omega_mom[3], 2) - pow(omega_mom[0], 2) - pow(omega_mom[1], 2) - pow(omega_mom[2], 2));

        value[i] = inv_mass_omega - mOmega;
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (abs(value[0]) < abs(value[1]))
      value_min = value[0];
    else
      value_min = value[1];

    return value_min;
  }

  Double_t minv_pi01_consv(Double_t *x, Double_t *p)
  {
    Int_t
        clusterind[4],
        clusterindpi0[2][2];

    Float_t
        gamma_mom[4][8],
        gamma_mom_pi0[2][2][4],
        pi0_mom[2][4],
        omega_mom[4],
        neu_vtx[4],
        trk[2][4],
        inv_mass_pi0;

    Double_t value[2], value_min, T0;

    T0 = p[30] * Trf;

    Reconstructor R;
    Solution S;

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3] + T0,
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }

    S = R.MySolve(selected);

    for (Int_t i = 0; i < 2; i++)
    {

      if (!S.error[i])
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
        neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
        neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
        neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

        Pi0PhotonPair(clusterind, gamma_mom, clusterindpi0, gamma_mom_pi0, pi0_mom, false, trk, omega_mom);

        inv_mass_pi0 = sqrt(pow(pi0_mom[0][3], 2) - pow(pi0_mom[0][0], 2) - pow(pi0_mom[0][1], 2) - pow(pi0_mom[0][2], 2));

        value[i] = inv_mass_pi0 - mPi0;
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (abs(value[0]) < abs(value[1]))
      value_min = value[0];
    else
      value_min = value[1];

    return value_min;
  }

  Double_t minv_pi02_consv(Double_t *x, Double_t *p)
  {
    Int_t
        clusterind[4],
        clusterindpi0[2][2];

    Float_t
        gamma_mom[4][8],
        gamma_mom_pi0[2][2][4],
        pi0_mom[2][4],
        omega_mom[4],
        neu_vtx[4],
        trk[2][4],
        inv_mass_pi0;

    Double_t value[2], value_min, T0;

    T0 = p[30] * Trf;

    Reconstructor R;
    Solution S;

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3] + T0,
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }

    S = R.MySolve(selected);

    for (Int_t i = 0; i < 2; i++)
    {

      if (!S.error[i])
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
        neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
        neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
        neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

        Pi0PhotonPair(clusterind, gamma_mom, clusterindpi0, gamma_mom_pi0, pi0_mom, false, trk, omega_mom);

        inv_mass_pi0 = sqrt(pow(pi0_mom[1][3], 2) - pow(pi0_mom[1][0], 2) - pow(pi0_mom[1][1], 2) - pow(pi0_mom[1][2], 2));

        value[i] = inv_mass_pi0 - mPi0;
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (abs(value[0]) < abs(value[1]))
      value_min = value[0];
    else
      value_min = value[1];

    return value_min;
  }

  Double_t x_consv(Double_t *x, Double_t *p)
  {
    Float_t gamma_mom[4][8], kaon_mom[4], neu_vtx[4], bhabha_vtx[3] = {p[24], p[25], p[26]},
                                                      y_axis[3] = {0., p[21], 0.}, ip[3], dist[3], kaon_vel[3], tot_length, tot_vel;

    Double_t value[2], value_min, T0;

    T0 = p[30] * Trf;

    Reconstructor R;
    Solution S;

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3] + T0,
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }

    S = R.MySolve(selected);

    for (Int_t i = 0; i < 2; i++)
    {

      if (S.error[i] == 0)
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        value[i] = neu_vtx[0] - p[20];
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (abs(value[0]) < abs(value[1]))
      value_min = value[0];
    else
      value_min = value[1];

    return value_min;
  }

  Double_t y_consv(Double_t *x, Double_t *p)
  {
    Float_t gamma_mom[4][8], kaon_mom[4], neu_vtx[4], bhabha_vtx[3] = {p[24], p[25], p[26]},
                                                      y_axis[3] = {0., p[21], 0.}, ip[3], dist[3], kaon_vel[3];

    Double_t value[2], value_min, T0;

    T0 = p[30] * Trf;

    Reconstructor R;
    Solution S;

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3] + T0,
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }

    S = R.MySolve(selected);

    for (Int_t i = 0; i < 2; i++)
    {

      if (S.error[i] == 0)
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        value[i] = neu_vtx[1] - p[21];
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (abs(value[0]) < abs(value[1]))
      value_min = value[0];
    else
      value_min = value[1];

    return value_min;
  }

  Double_t z_consv(Double_t *x, Double_t *p)
  {
    Float_t gamma_mom[4][8], kaon_mom[4], neu_vtx[4], bhabha_vtx[3] = {p[24], p[25], p[26]},
                                                      y_axis[3] = {0., p[21], 0.}, ip[3], dist[3], kaon_vel[3];

    Double_t value[2], value_min, T0;

    T0 = p[30] * Trf;

    Reconstructor R;
    Solution S;

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3] + T0,
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }

    S = R.MySolve(selected);

    for (Int_t i = 0; i < 2; i++)
    {

      if (S.error[i] == 0)
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        value[i] = neu_vtx[2] - p[22];
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (abs(value[0]) < abs(value[1]))
      value_min = value[0];
    else
      value_min = value[1];

    return value_min;
  }

}