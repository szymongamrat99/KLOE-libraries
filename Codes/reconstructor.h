#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

struct Solution{
  /*
    1,2
    x,y,z,t
  */
  double sol[2][4];
  // indicates error computing solution
  bool error[2];
};

class Reconstructor{
  
 public:
  void SetClu(int i, double x, double y, double z, double t, double E);
  Solution MySolve(int * selected);
  Solution KleusbergSolve(int * selected);
  Solution LeastSquaresSolve(double * x0);
  Solution MinuitSolve(double *x0);
  
  double ResidualErr(int i, double const * x);
  double ResidualErrTot(double const * x);
  double CombEnergy(int * selected);
  double TotalEnergy()const;
  double GetInvMasses(double * sol, int * comb, double * Minvgg)const;
  double GetInvMassDiscrepancy(double * sol, int * comb)const;
  void GetKmomentum(const double * sol, float * p)const;
  
 private:
  /*
    x,y,z,t
   */
  double _clu[6][4];
  double _ene[6];
};

#endif
