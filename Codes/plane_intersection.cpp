#include "plane_intersection.h"
#include "TMath.h"
#include "../const.h"

void plane_intersection(Float_t pointplane[3], Float_t vectorplane1[3], Float_t vectorplane2[3],
                        Float_t pointline[3], Float_t vectorline[3],
                        Float_t ip_plane[3])
{
    Float_t n[3], d, l_n, numerator;

    n[0] = vectorplane1[1]*vectorplane2[2] - vectorplane1[2]*vectorplane2[1];
    n[1] = vectorplane1[2]*vectorplane2[0] - vectorplane1[0]*vectorplane2[2];
    n[2] = vectorplane1[0]*vectorplane2[1] - vectorplane1[1]*vectorplane2[0];
    
    l_n = vectorline[0]*n[0] + vectorline[1]*n[1] + vectorline[2]*n[2];

    numerator = (pointplane[0] - pointline[0])*n[0] +
                (pointplane[1] - pointline[1])*n[1] +
                (pointplane[2] - pointline[2])*n[2];

    d = numerator/l_n;

    ip_plane[0] = pointline[0] + vectorline[0]*d;
    ip_plane[1] = pointline[1] + vectorline[1]*d;
    ip_plane[2] = pointline[2] + vectorline[2]*d;
}