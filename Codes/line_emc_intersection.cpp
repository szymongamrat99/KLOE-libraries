// Calculation of intersection of line including 
// neutral vtx and directed with gamma momentum
// With the EMC

#include "line_emc_intersection.h"

Int_t cluster_finder(Float_t *point, Float_t *vector, Float_t *cluster)
{
  TVector3 neu_vtx(point), pgamma(vector), z_axis(0., 0., 1.), neg_vec(0., 0., zmin), pos_vec(0., 0., zmax);
  TVector3 neu_vtx_T(point), pgamma_T(vector); // Transverse components

  neu_vtx_T(2) = 0.;
  pgamma_T(2) = 0.;

  Float_t dot_pr = 0., cross_pr = 0., m2_vec = 0., m2_pnt = 0., sq_root = 0.;
  Float_t numerator[2] = {0., 0.}, numerator_plane = 0., denominator = 0., t[2] = {0., 0.}, t_plane = 0.;

  // First determine the parameter of a line, for which it intersects with the barrel
  m2_vec = pgamma_T.Mag2();
  denominator = m2_vec;
  //
  //
  dot_pr = pgamma_T.Dot(neu_vtx_T);
  cross_pr = (pgamma_T.Cross(neu_vtx_T)).Mag2();

  sq_root = sqrt(m2_vec * pow(RT, 2) - cross_pr);

  numerator[0] = -dot_pr + sq_root;
  numerator[1] = -dot_pr - sq_root;

  t[0] = numerator[0] / denominator;
  t[1] = numerator[1] / denominator;

  // Check, which point is parallel to the momentum
  TVector3 barrel_points[2] = {neu_vtx + pgamma * t[0], neu_vtx + pgamma * t[1]};
  Float_t angle[2] = {0., 0.};

  angle[0] = (barrel_points[0] - neu_vtx).Angle(pgamma);
  angle[1] = (barrel_points[1] - neu_vtx).Angle(pgamma);

  TVector3 inter_point;

  if (angle[0] < 0.01 && angle[1] > 3.)
    inter_point.SetXYZ(barrel_points[0](0), barrel_points[0](1), barrel_points[0](2));
  else if (angle[1] < 0.01 && angle[0] > 3.)
    inter_point.SetXYZ(barrel_points[1](0), barrel_points[1](1), barrel_points[1](2));

  cluster[0] = inter_point(0);
  cluster[1] = inter_point(1);
  cluster[2] = inter_point(2);

  // Check if |z| < 165 cm, if not find intersection with the plane

  if (inter_point(2) > 165)
  {
    numerator_plane = (pos_vec - neu_vtx).Dot(z_axis);
    denominator = pgamma.Dot(z_axis);

    t_plane = numerator_plane / denominator;

    inter_point = neu_vtx + pgamma * t_plane;

    cluster[0] = inter_point(0);
    cluster[1] = inter_point(1);
    cluster[2] = inter_point(2);

    return 1;
  }
  else if (inter_point(2) < -165)
  {
    numerator_plane = (neg_vec - neu_vtx).Dot(z_axis);
    denominator = pgamma.Dot(z_axis);

    t_plane = numerator_plane / denominator;

    inter_point = neu_vtx + pgamma * t_plane;

    cluster[0] = inter_point(0);
    cluster[1] = inter_point(1);
    cluster[2] = inter_point(2);

    return 1;
  }
  else
    return 0;
}